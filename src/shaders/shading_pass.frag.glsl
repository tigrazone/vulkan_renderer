//  Copyright (C) 2021, Christoph Peters, Karlsruhe Institute of Technology
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <https://www.gnu.org/licenses/>.

#version 460
#extension GL_GOOGLE_include_directive : enable
#extension GL_EXT_samplerless_texture_functions : enable
#extension GL_EXT_nonuniform_qualifier : enable
#extension GL_EXT_control_flow_attributes : enable

#if RTX_ON
#extension GL_EXT_ray_query : enable
#endif

#include "shared_constants.glsl"

#include "noise_utility.glsl"
#include "brdfs.glsl"
#include "mesh_quantization.glsl"
#include "polygon_sampling.glsl"
#include "polygon_clipping.glsl"
#include "srgb_utility.glsl"

//! Bindings for mesh geometry (see mesh_t in the C code)
layout(binding = 1) uniform utextureBuffer g_quantized_vertex_positions;
layout(binding = 2) uniform textureBuffer g_packed_normals_and_tex_coords;
layout(binding = 3) uniform utextureBuffer g_material_indices;

//! The texture with primitive indices per pixel produced by the visibility pass
layout(binding = 4, input_attachment_index = 0) uniform usubpassInput g_visibility_buffer;

//! Textures (base color, specular, normal consecutively) for each material
//layout(binding = 5) uniform sampler2D g_material_textures[3 * MATERIAL_COUNT];
layout(binding = 5) uniform sampler2D g_material_textures[];

//! Textures for each polygonal light. These can be plane space textures, light
//! probes or IES profiles
//layout(binding = 8) uniform sampler2D g_light_textures[LIGHT_TEXTURE_COUNT];
layout(binding = 8) uniform sampler2D g_light_textures[];

//! The pixel index with origin in the upper left corner
layout(origin_upper_left) in vec4 gl_FragCoord;
//! Color written to the swapchain image
layout(location = 0) out vec4 g_out_color;

#if RTX_ON
//! The top-level acceleration structure that contains all shadow-casting
//! geometry
layout(binding = 9, set = 0) uniform accelerationStructureEXT g_top_level_acceleration_structure;
#else
//! Software raytracing data
#endif

/*! If shadow rays are enabled, this function traces a shadow ray towards the
	given polygonal light and updates visibility accordingly. If visibility is
	false already, no ray is traced. The ray direction must be normalized.*/
void get_polygon_visibility(inout bool visibility, vec3 sampled_dir, vec3 shading_position, polygonal_light_t polygonal_light)
{
	if ((g_rtx_bits & rtx_bits_TRACE_SHADOW_RAYS) > 0)
	{
#if RTX_ON
		if (visibility)
		{
			float max_t = -dot(vec4(shading_position, 1.0f), polygonal_light.plane) / dot(sampled_dir, polygonal_light.plane.xyz);
			float min_t = 1.0e-3f;
			// Perform a ray query and wait for it to finish. One call to
			// rayQueryProceedEXT() should be enough because of
			// gl_RayFlagsTerminateOnFirstHitEXT.
			rayQueryEXT ray_query;
			rayQueryInitializeEXT(ray_query, g_top_level_acceleration_structure,
								  gl_RayFlagsTerminateOnFirstHitEXT | gl_RayFlagsOpaqueEXT | gl_RayFlagsSkipClosestHitShaderEXT,
								  0xFF, shading_position, min_t, sampled_dir, max_t);
			rayQueryProceedEXT(ray_query);
			// Update the visibility
			bool occluder_hit = (rayQueryGetIntersectionTypeEXT(ray_query, true) != gl_RayQueryCommittedIntersectionNoneEXT);
			visibility = !occluder_hit;
		}
#else
		// software raytracing
#endif
	}
}

/*! Determines the radiance received from the given direction due to the given
	polygonal light (ignoring visibility).
	\param sampled_dir The normalized direction from the shading point to the
		light source in world space. It must be chosen such that the ray
		actually intersects the polygon (possibly behind an occluder or below
		the horizon).
	\param shading_position The location of the shading point.
	\param polygonal_light The light source for which the incoming radiance is
		evaluated.
	\return Received radiance.*/
vec3 get_polygon_radiance(vec3 sampled_dir, vec3 shading_position, polygonal_light_t polygonal_light)
{
	vec3 radiance = polygonal_light.surface_radiance;
	uint technique = polygonal_light.texturing_technique;
	if (technique != polygon_texturing_none)
	{
		vec2 tex_coord;
		if (technique == polygon_texturing_area)
		{
			// Intersect the ray with the plane of the light source
			float intersection_t = -dot(vec4(shading_position, 1.0f), polygonal_light.plane) / dot(sampled_dir, polygonal_light.plane.xyz);
			vec3 intersection = shading_position + intersection_t * sampled_dir;
			// Transform to plane space
			intersection -= polygonal_light.translation;
			tex_coord = (transpose(polygonal_light.rotation) * intersection).xy;
			tex_coord *= vec2(polygonal_light.inv_scaling_x, polygonal_light.inv_scaling_y);
		}
		else
		{
			vec3 lookup_dir;
			if (technique == polygon_texturing_ies_profile)
			{
				// For IES profiles, we transform to plane space
				lookup_dir = transpose(polygonal_light.rotation) * sampled_dir;
				// IES profiles already include this cosine term, so we divide
				// it out
				radiance /= abs(lookup_dir.z);
			}
			else
				// This code is designed to be compatible with the coordinate
				// system used for HDRI Haven light probes
				lookup_dir = vec3(-sampled_dir.x, sampled_dir.y, sampled_dir.z);
			// Now we compute spherical coordinates
			tex_coord.x = atan(lookup_dir.y, lookup_dir.x) * M_HALF_INV_PI;
			tex_coord.y = acos(lookup_dir.z) * M_INV_PI;
		}
		radiance *= textureLod(g_light_textures[polygonal_light.texture_index], tex_coord, 0.0f).rgb;
	}
	return radiance;
}

/*! Determines the radiance received from the given direction due to the given
	polygonal light and multiplies it by the BRDF for this direction. If
	necessary, this function traces a shadow ray to determine visibility.
	\param lambert The dot product of the normal vector and the given
		direction, in case you have use for it outside this function.
	\param sampled_dir The normalized direction from the shading point to the
		light source in world space. It must be chosen such that the ray
		actually intersects the polygon (possibly behind an occluder or below
		the horizon).
	\param shading_data Shading data for the shading point.
	\param polygonal_light The light source for which the incoming radiance is
		evaluated.
	\param diffuse Whether the diffuse BRDF component is evaluated.
	\param specular Whether the specular BRDF component is evaluated.
	\return BRDF times incoming radiance times visibility.*/

//! Like get_polygon_radiance_visibility_brdf_product() but always evaluates
//! both BRDF components and also outputs the visibility term explicitly.
vec3 get_polygon_radiance_visibility_brdf_product(out bool out_visibility, vec3 sampled_dir, shading_data_t shading_data, polygonal_light_t polygonal_light)
{
	out_visibility = (dot(shading_data.normal, sampled_dir) > 0.0f);
	if (out_visibility)
	{
		get_polygon_visibility(out_visibility, sampled_dir, shading_data.position, polygonal_light);
		if (!out_visibility)
			return vec3(0.0f);
		vec3 radiance = get_polygon_radiance(sampled_dir, shading_data.position, polygonal_light);
		return radiance * evaluate_brdf(shading_data, sampled_dir);
	}
	else
		return vec3(0.0f);
}

/*! Returns the MIS estimator for the given sample using the currently enabled
	MIS heuristic. It supports our weighted balance heuristic and optimal MIS.
	\param visibility true iff the given sample is occluded.
	\param integrand The value of the integrand with respect to solid angle
		measure at the sampled location.
	\param sampled_density, other_density See get_mis_weight_over_density().
	\param sampled_weight, other_weight Estimates of unshadowed shading for the
		respective BRDF components. For all techniques except optimal MIS, it
		is legal to introduce an arbitrary but identical constant factor for
		both of them.
	\param visibility_estimate Some estimate of how much the shading point is
		shadowed on average (0 means fully shadowed). It does not have to be
		accurate, it just blends between two MIS heuristics for optimal MIS.
	\return An unbiased multiple importance sampling estimator. Note that it
		may be negative when optimal MIS is used.*/
vec3 get_mis_estimate(bool visibility, vec3 integrand, vec3 sampled_weight, float sampled_density, vec3 other_weight, float other_density, float visibility_estimate)
{
	float balance_weight_over_density = 1.0f / (sampled_density + other_density);
	vec3 weighted_sum = sampled_weight * sampled_density + other_weight * other_density;
	vec3 weighted_weight_over_density = sampled_weight / weighted_sum;
	vec3 mixed_weight_over_density = vec3(fma(-visibility_estimate, balance_weight_over_density, balance_weight_over_density));
	mixed_weight_over_density = fma(vec3(visibility_estimate), weighted_weight_over_density, vec3(mixed_weight_over_density));
	// For visible samples, we use the actual integrand
	vec3 visible_estimate = mixed_weight_over_density * integrand;
	return visible_estimate;
}

/*! Takes samples from the given polygonal light to compute shading. The number
	of samples and sampling techniques are determined by defines.
	\return The color that arose from shading.*/
vec3 evaluate_polygonal_light_shading(shading_data_t shading_data, ltc_coefficients_t ltc, polygonal_light_t polygonal_light, inout noise_accessor_t accessor)
{
	vec3 result = vec3(0.0f);

	// If the shading point is on the wrong side of the polygon, we get a
	// correct winding by flipping the orientation of the shading space
	float side = dot(vec4(shading_data.position, 1.0f), polygonal_light.plane);
	if (side < 0.0f)
	{
		[[unroll]] for (uint i = 0; i != 4; ++i)
		{
			ltc.world_to_shading_space[i][1] = -ltc.world_to_shading_space[i][1];
			ltc.world_to_cosine_space[i][1] = -ltc.world_to_cosine_space[i][1];
		}
	}

	// Combined diffuse and specular strategies

	// Instruction cache misses are a concern. Thus, we strive to keep the code
	// small by preparing the diffuse (i==0) and specular (i==1) sampling
	// strategies in the same loop.
	projected_solid_angle_polygon_t polygon_diffuse;
	projected_solid_angle_polygon_t polygon_specular;

	vec3 vertices_local_space[MAX_POLYGON_VERTEX_COUNT];

	[[dont_unroll]] for (uint i = 0; i != 2; ++i)
	{
		// Local space is either shading space (for the diffuse technique) or
		// cosine space (for the specular technique)
		mat4x3 world_to_local_space = (i == 0) ? ltc.world_to_shading_space : ltc.world_to_cosine_space;
		if (i > 0)
			// We put this object in the wrong place at first to avoid move
			// instructions
			polygon_diffuse = polygon_specular;
		// Transform to local space
		[[unroll]] for (uint j = 0; j != MAX_POLYGONAL_LIGHT_VERTEX_COUNT; ++j)
			vertices_local_space[j] = world_to_local_space * vec4(polygonal_light.vertices_world_space[j], 1.0f);
		// Clip
		uint clipped_vertex_count = clip_polygon(polygonal_light.vertex_count, vertices_local_space);
		if (clipped_vertex_count == 0 && i == 0)
			// The polygon is completely below the horizon
			return vec3(0.0f);
		else if (clipped_vertex_count == 0)
		{
			// The linearly transformed cosine is zero on the polygon
			polygon_specular.projected_solid_angle = 0.0f;
			break;
		}
		// Prepare sampling
		polygon_specular = prepare_projected_solid_angle_polygon_sampling(clipped_vertex_count, vertices_local_space);
	}
	// Even when something remains after clipping, the projected solid angle
	// may still underflow
	if (polygon_diffuse.projected_solid_angle == 0.0f)
		return vec3(0.0f);
	// Compute the importance of the specular sampling technique using an
	// LTC-based estimate of unshadowed shading
	float specular_albedo = ltc.albedo;
	float specular_weight = specular_albedo * polygon_specular.projected_solid_angle;

	// Compute the importance of the diffuse sampling technique using the
	// diffuse albedo and the projected solid angle. Zero albedo is forbidden
	// because we need a non-zero weight for diffuse samples in parts where the
	// LTC is zero but the specular BRDF is not. Thus, we clamp.
	vec3 diffuse_albedo = max(shading_data.diffuse_albedo, vec3(0.01f));
	vec3 diffuse_weight = diffuse_albedo * polygon_diffuse.projected_solid_angle;
	uint technique_count = (polygon_specular.projected_solid_angle > 0.0f) ? 2 : 1;
	float rcp_diffuse_projected_solid_angle = 1.0f / polygon_diffuse.projected_solid_angle;
	float rcp_specular_projected_solid_angle = 1.0f / polygon_specular.projected_solid_angle;
	vec3 specular_weight_rgb = vec3(specular_weight);
	// For optimal MIS, constant factors in the diffuse and specular weight
	// matter

	// Take the requested number of samples with both techniques
	for (uint i = 0; i != SAMPLE_COUNT; ++i) {
		// Take the samples
		vec3 dir_shading_space_diffuse = sample_projected_solid_angle_polygon(polygon_diffuse, get_noise_2(accessor));
		vec3 dir_shading_space_specular;
		if (polygon_specular.projected_solid_angle > 0.0f)
		{
			dir_shading_space_specular = sample_projected_solid_angle_polygon(polygon_specular, get_noise_2(accessor));
			dir_shading_space_specular = normalize(ltc.cosine_to_shading_space * dir_shading_space_specular);
		}
		[[dont_unroll]] for (uint j = 0; j != technique_count; ++j)
		{
			vec3 dir_shading_space = (j == 0) ? dir_shading_space_diffuse : dir_shading_space_specular;
			if (dir_shading_space.z <= 0.0f)
				continue;
			// Compute the densities for the sample with respect to both
			// sampling techniques (w.r.t. solid angle measure)
			float diffuse_density = dir_shading_space.z * rcp_diffuse_projected_solid_angle;
			float specular_density = evaluate_ltc_density(ltc, dir_shading_space, rcp_specular_projected_solid_angle);
			// Evaluate radiance and BRDF and the integrand as a whole
			bool visibility;
			vec3 integrand = dir_shading_space.z * get_polygon_radiance_visibility_brdf_product(visibility, (transpose(ltc.world_to_shading_space) * dir_shading_space).xyz, shading_data, polygonal_light);
			// Use the appropriate MIS heuristic to turn the sample into a
			// splat and accummulate
			if (j == 0 && polygon_specular.projected_solid_angle <= 0.0f)
			{
				// We only have one sampling technique, so no MIS is needed
				if (visibility)
					result += integrand / diffuse_density;
			}
			else if (j == 0)
				result += get_mis_estimate(visibility, integrand, diffuse_weight, diffuse_density, specular_weight_rgb, specular_density, g_mis_visibility_estimate);
			else
				result += get_mis_estimate(visibility, integrand, specular_weight_rgb, specular_density, diffuse_weight, diffuse_density, g_mis_visibility_estimate);
		}
	}

	return result / SAMPLE_COUNT;
}

/*! Based on the knowledge that the given primitive is visible on the given
	pixel, this function recovers complete shading data for this pixel.
	\param pixel Coordinates of the pixel inside the viewport in pixels.
	\param primitive_index Index into g_material_indices, etc.
	\param ray_direction Direction of a ray through that pixel in world space.
		Does not need to be normalized.
	\note This implementation assumes a perspective projection */
shading_data_t get_shading_data(ivec2 pixel, int primitive_index, vec3 ray_direction)
{
	shading_data_t result;
	// Load position, normal and texture coordinates for each triangle vertex
	vec3 positions[3], normals[3];
	vec2 tex_coords[3];
	int vertex_index = primitive_index + primitive_index + primitive_index;
	[[unroll]] for (int i = 0; i != 3; ++i)
	{
		uvec2 quantized_position = texelFetch(g_quantized_vertex_positions, vertex_index).rg;
		positions[i] = decode_position_64_bit(quantized_position, g_mesh_dequantization_factor, g_mesh_dequantization_summand);
		vec4 normal_and_tex_coords = texelFetch(g_packed_normals_and_tex_coords, vertex_index);
		normals[i] = decode_normal_32_bit(normal_and_tex_coords.xy);
		tex_coords[i] = fma(normal_and_tex_coords.zw, vec2(8.0f, -8.0f), vec2(0.0f, 1.0f));
		vertex_index++;
	}
	// Construct the view ray for the pixel at hand (the ray direction is not
	// normalized)
	vec3 ray_origin = g_camera_position_world_space;
	// Perform ray triangle intersection to figure out barycentrics within the
	// triangle
	vec3 barycentrics;
	vec3 edges[2] = {
		positions[1] - positions[0],
		positions[2] - positions[0]};
	vec3 ray_cross_edge_1 = cross(ray_direction, edges[1]);
	float rcp_det_edges_direction = 1.0f / dot(edges[0], ray_cross_edge_1);
	vec3 ray_to_0 = ray_origin - positions[0];
	float det_0_dir_edge_1 = dot(ray_to_0, ray_cross_edge_1);
	barycentrics.y = rcp_det_edges_direction * det_0_dir_edge_1;
	vec3 edge_0_cross_0 = cross(edges[0], ray_to_0);
	float det_dir_edge_0_0 = dot(ray_direction, edge_0_cross_0);
	barycentrics.z = -rcp_det_edges_direction * det_dir_edge_0_0;
	barycentrics.x = 1.0f - (barycentrics.y + barycentrics.z);
	// Compute screen space derivatives for the barycentrics
	vec3 barycentrics_derivs[2];
	[[unroll]] for (uint i = 0; i != 2; ++i)
	{
		vec3 ray_direction_deriv = g_pixel_to_ray_direction_world_space[i];
		vec3 ray_cross_edge_1_deriv = cross(ray_direction_deriv, edges[1]);
		float rcp_det_edges_direction_deriv = -dot(edges[0], ray_cross_edge_1_deriv) * rcp_det_edges_direction * rcp_det_edges_direction;
		float det_0_dir_edge_1_deriv = dot(ray_to_0, ray_cross_edge_1_deriv);
		barycentrics_derivs[i].y = rcp_det_edges_direction_deriv * det_0_dir_edge_1 + rcp_det_edges_direction * det_0_dir_edge_1_deriv;
		float det_dir_edge_0_0_deriv = dot(ray_direction_deriv, edge_0_cross_0);
		barycentrics_derivs[i].z = -rcp_det_edges_direction_deriv * det_dir_edge_0_0 - rcp_det_edges_direction * det_dir_edge_0_0_deriv;
		barycentrics_derivs[i].x = -(barycentrics_derivs[i].y + barycentrics_derivs[i].z);
	}
	// Interpolate vertex attributes across the triangle
	result.position = fma(vec3(barycentrics[0]), positions[0], fma(vec3(barycentrics[1]), positions[1], barycentrics[2] * positions[2]));
	vec3 interpolated_normal = normalize(fma(vec3(barycentrics[0]), normals[0], fma(vec3(barycentrics[1]), normals[1], barycentrics[2] * normals[2])));
	vec2 tex_coord = fma(vec2(barycentrics[0]), tex_coords[0], fma(vec2(barycentrics[1]), tex_coords[1], barycentrics[2] * tex_coords[2]));
	// Compute screen space texture coordinate derivatives for filtering
	vec2 tex_coord_derivs[2] = {vec2(0.0f), vec2(0.0f)};
	[[unroll]] for (uint i = 0; i != 2; ++i)
		[[unroll]] for (uint j = 0; j != 3; ++j)
			tex_coord_derivs[i] += barycentrics_derivs[i][j] * tex_coords[j];
	// Read all three textures
	uint material_index = texelFetch(g_material_indices, primitive_index).r;
	uint material_index3 = material_index + material_index + material_index;
	vec3 base_color = textureGrad(g_material_textures[nonuniformEXT(material_index3)], tex_coord, tex_coord_derivs[0], tex_coord_derivs[1]).rgb;
	vec3 specular_data = textureGrad(g_material_textures[nonuniformEXT(material_index3 + 1)], tex_coord, tex_coord_derivs[0], tex_coord_derivs[1]).rgb;
	vec3 normal_tangent_space;
	normal_tangent_space.xy = textureGrad(g_material_textures[nonuniformEXT(material_index3 + 2)], tex_coord, tex_coord_derivs[0], tex_coord_derivs[1]).rg;
	normal_tangent_space.xy = fma(normal_tangent_space.xy, vec2(2.0f), vec2(-1.0f));
	normal_tangent_space.z = sqrt(max(0.0f, fma(-normal_tangent_space.x, normal_tangent_space.x, fma(-normal_tangent_space.y, normal_tangent_space.y, 1.0f))));
	// Prepare BRDF parameters (i.e. immitate Falcor to be compatible with its
	// assets, which in turn immitates the Unreal engine). The Fresnel F0 value
	// for surfaces with zero metalicity is set to 0.02, not 0.04 as in Falcor
	// because this way colors throughout the scenes are a little less
	// desaturated.
	float metalicity = specular_data.b;
	result.diffuse_albedo = fma(base_color, -vec3(metalicity), base_color);
	result.fresnel_0 = mix(vec3(0.02f), base_color, metalicity);
	float linear_roughness = specular_data.g;
	result.roughness = linear_roughness * linear_roughness;
	result.roughness = clamp(result.roughness * g_roughness_factor, 0.0064f, 1.0f);
	// Transform the normal vector to world space
	vec2 tex_coord_edges[2] = {
		tex_coords[1] - tex_coords[0],
		tex_coords[2] - tex_coords[0]};
	vec3 normal_cross_edge_0 = cross(interpolated_normal, edges[0]);
	vec3 edge1_cross_normal = cross(edges[1], interpolated_normal);
	vec3 tangent = edge1_cross_normal * tex_coord_edges[0].x + normal_cross_edge_0 * tex_coord_edges[1].x;
	vec3 bitangent = edge1_cross_normal * tex_coord_edges[0].y + normal_cross_edge_0 * tex_coord_edges[1].y;
	float mean_tangent_length = sqrt(0.5f * (dot(tangent, tangent) + dot(bitangent, bitangent)));
	mat3 tangent_to_world_space = mat3(tangent, bitangent, interpolated_normal);
	normal_tangent_space.z *= max(1.0e-10f, mean_tangent_length);
	result.normal = normalize(tangent_to_world_space * normal_tangent_space);
	// Perform local shading normal adaptation to avoid that the view direction
	// is below the horizon. Inspired by Keller et al., Section A.3, but
	// different since the method of Keller et al. often led to normal vectors
	// "running off to the side." We simply clip the shading normal into the
	// hemisphere of the outgoing direction.
	// https://arxiv.org/abs/1705.01263
	result.outgoing = normalize(g_camera_position_world_space - result.position);
	float normal_offset = max(0.0f, 1.0e-3f - dot(result.normal, result.outgoing));
	result.normal = fma(vec3(normal_offset), result.outgoing, result.normal);
	result.normal = normalize(result.normal);
	result.lambert_outgoing = dot(result.normal, result.outgoing);
	return result;
}

void main()
{
	// Obtain an integer pixel index
	ivec2 pixel = ivec2(gl_FragCoord.xy);
	// Get the primitive index from the visibility buffer
	uint primitive_index = subpassLoad(g_visibility_buffer).r;
	// Set the background color
	vec3 final_color = vec3(0.0f);
	// Figure out the ray to the first visible surface
	vec3 view_ray_direction = g_pixel_to_ray_direction_world_space * vec3(pixel, 1.0f);
	vec4 view_ray_end;
	shading_data_t shading_data;
	if (primitive_index == 0xFFFFFFFF)
		view_ray_end = vec4(view_ray_direction, 0.0f);
	else
	{
		// Prepare shading data for the visible surface point
		shading_data = get_shading_data(pixel, int(primitive_index), view_ray_direction);
		view_ray_end = vec4(shading_data.position, 1.0f);
	}

	if(SHOW_POLYGONAL_LIGHTS == 1)
	{
		// Display light sources
		view_ray_direction = normalize(view_ray_direction);
		for (uint i = 0; i != POLYGONAL_LIGHT_COUNT; ++i)
			if (polygonal_light_ray_intersection(g_polygonal_lights[i], g_camera_position_world_space, view_ray_end))
				final_color += get_polygon_radiance(view_ray_direction, g_camera_position_world_space, g_polygonal_lights[i]);
	}

	// We only need to shade anything if there is a primitive to shade
	if (primitive_index != 0xFFFFFFFF)
	{
		// Get ready to use linearly transformed cosines
		float fresnel_luminance = dot(shading_data.fresnel_0, luminance_weights);
		ltc_coefficients_t ltc = get_ltc_coefficients(fresnel_luminance, shading_data.roughness, shading_data.position, shading_data.normal, shading_data.outgoing, g_ltc_constants);
		// Prepare noise for all sampling decisions
		noise_accessor_t noise_accessor = get_noise_accessor(pixel, g_noise_resolution_mask, g_noise_texture_index_mask, g_noise_random_numbers);
		// Shade with all polygonal lights
		for (uint i = 0; i != POLYGONAL_LIGHT_COUNT; ++i) {
			final_color += evaluate_polygonal_light_shading(shading_data, ltc, g_polygonal_lights[i], noise_accessor);
		}
	}
	// If there are NaNs or INFs, we want to know. Make them pink.
	if (isnan(final_color.r) || isnan(final_color.g) || isnan(final_color.b) || isinf(final_color.r) || isinf(final_color.g) || isinf(final_color.b))
		final_color = vec3(1.0f, 0.0f, 0.8f) / g_exposure_factor;
	// Output the result of shading
	g_out_color = vec4(final_color * g_exposure_factor, 1.0f);
	// Here is how we support HDR screenshots: Always rendering to an
	// intermediate HDR render target would be wasting resources, since we do
	// not take screenshots each frame. Instead, a HDR screenshot consists of
	// two LDR screenshots holding different bits of halfs.
	if (g_frame_bits > 0)
	{
		uint mask = (g_frame_bits == 1) ? 0xFF : 0xFF00;
		uint shift = (g_frame_bits == 1) ? 0 : 8;
		uvec2 half_bits = uvec2(packHalf2x16(g_out_color.rg), packHalf2x16(g_out_color.ba));
		g_out_color = vec4(
			//((half_bits[0] & mask) >> shift) * (1.0f / 255.0f),
			((half_bits[0] & mask) >> shift) * INV_255,
			((((half_bits[0] & 0xFFFF0000) >> 16) & mask) >> shift) * INV_255,
			((half_bits[1] & mask) >> shift) * INV_255,
			1.0f);
		// We just want to write bits to the render target, not colors. If the
		// graphics pipeline does linear to sRGB conversion for us, we do sRGB
		// to linear conversion here to counter that.
		if(OUTPUT_LINEAR_RGB> 0) {
			g_out_color.rgb = convert_srgb_to_linear_rgb(g_out_color.rgb);
		}
	}
	// if g_frame_bits == 0, we output linear RGB or sRGB as requested
	else
		if(OUTPUT_LINEAR_RGB == 0) {
			g_out_color.rgb = convert_linear_rgb_to_srgb(g_out_color.rgb);
		}
}