/****************************************************
 * 2016-2021 Eric Bachmann and Mike Zmuda
 * All Rights Reserved.
 * NOTICE:
 * Dissemination of this information or reproduction
 * of this material is prohibited unless prior written
 * permission is granted.
 ****************************************************/

#include "fragmentops.h"

FogParams FragmentOps::fogParams;
bool FragmentOps::performDepthTest = true;
bool FragmentOps::readonlyDepthBuffer = false;
bool FragmentOps::readonlyColorBuffer = false;

/**
 * @fn	double FogParams::fogFactor(const dvec3 &fragPos, const dvec3 &eyePos) const
 * @brief	Computes fog factor - f.
 * @param	fragPos	The fragment position.
 * @param	eyePos 	The eye position.
 * @return	The fog factor - f.
 */

double FogParams::fogFactor(const dvec3 &fragPos, const dvec3 &eyePos) const {
	/* CSE 386 - todo  */
	double fogFactor = 1.0;
	return fogFactor;
}

/**
 * @fn	color FragmentOps::applyLighting(const Fragment &fragment, 
 *										const dvec3 &eyePositionInWorldCoords,
 *										const vector<LightSourcePtr> &lights,
 *										const dmat4 &viewingMatrix)
 * @brief	Applies the lighting to a fragment
 * @param	fragment					The fragment.
 * @param	eyePositionInWorldCoords	The eye position in world coordinates.
 * @param	lights						The vector of lights in the scene.
 * @param               eyeFrame        The camera's frame.
 *
 * @return	The color of the fragment after applying lighting equations.
 */

color FragmentOps::applyLighting(const Fragment &fragment, const dvec3 &eyePositionInWorldCoords,
										const vector<LightSourcePtr> &lights,
										const Frame &eyeFrame) {
	/* CSE 386 - todo  */
	return fragment.material.diffuse;
}

/**
 * @fn	color FragmentOps::applyFog(const color &destColor, const dvec3 &eyePos, const dvec3 &fragPos)
 * @brief	Applies fog to a fragment.
 * @param	destColor	Destination color.
 * @param	eyePos   	Eye position.
 * @param	fragPos  	Fragment position.
 * @return	The color after applying the fog.
 */

color FragmentOps::applyFog(const color &destColor,
							const dvec3 &eyePos, const dvec3 &fragPos) {
	/* CSE 386 - todo  */
	return destColor;
}

/**
 * @fn	color FragmentOps::applyBlending(double alpha, const color &srcColor, const color &destColor)
 * @brief	Applies blending to a fragment.
 * @param	alpha	 	Alpha value.
 * @param	srcColor 	Source color.
 * @param	destColor	Destination color.
 * @return	The blended color.
 */

color FragmentOps::applyBlending(double alpha, const color &srcColor, const color &destColor) {
	/* CSE 386 - todo  */
	return srcColor;
}

/**
 * @fn	void FragmentOps::processFragment(FrameBuffer &frameBuffer, 
 *											const dvec3 &eyePositionInWorldCoords,
 *											const vector<LightSourcePtr> lights, 
 *											const Fragment &fragment,
 *											const dmat4 &viewingMatrix)
 * @brief	Process the fragment, leaving the results in the framebuffer.
 * @param [in,out]	frameBuffer	                The frame buffer
 * @param 		  	eyePositionInWorldCoords	The eye position in world coordinates.
 * @param 		  	lights						Vector of lights in scene.
 * @param 		  	fragment					Fragment to be processed.
 * @param           eyeFrame                    The camera's frame.
 */

void FragmentOps::processFragment(FrameBuffer& frameBuffer, const dvec3& eyePositionInWorldCoords,
	const vector<LightSourcePtr> lights,
	const Fragment& fragment,
	const Frame& eyeFrame) {
	const dvec3& eyePos = eyePositionInWorldCoords;

	double Z = fragment.windowPos.z;
	int X = (int)fragment.windowPos.x;
	int Y = (int)fragment.windowPos.y;
	DEBUG_PIXEL = (X == xDebug && Y == yDebug);

	/* CSE 386 - todo */
	color C = fragment.material.diffuse;
	frameBuffer.setColor(X, Y, C);
	frameBuffer.setDepth(X, Y, Z);
}