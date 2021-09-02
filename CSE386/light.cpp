/****************************************************
 * 2016-2021 Eric Bachmann and Mike Zmuda
 * All Rights Reserved.
 * NOTICE:
 * Dissemination of this information or reproduction
 * of this material is prohibited unless prior written
 * permission is granted.
 ****************************************************/

#include "light.h"
#include "io.h"

/**
 * @fn	color ambientColor(const color &mat, const color &light)
 * @brief	Computes the ambient color produced by a single light at a single point.
 * @param	mat  	Ambient material property.
 * @param	lightAmbient	Light's ambient color.
 * @return	Ambient color.
  */

color ambientColor(const color &mat, const color &lightAmbient) {
	/* CSE 386 - todo  */
	return mat;
}

/**
 * @fn	color diffuseColor(const color &mat, const color &light, const dvec3 &l, const dvec3 &n)
 * @brief	Computes diffuse color produce by a single light at a single point.
 * @param	mat		 	Material.
 * @param	lightDiffuse	 	The light color.
 * @param	l		 	Light vector.
 * @param	n		 	Normal vector.
 * @return	Diffuse color.
 */

color diffuseColor(const color &mat, const color &lightDiffuse,
					const dvec3 &l, const dvec3 &n) {
	/* CSE 386 - todo  */
	return mat;
}

/**
 * @fn	color specularColor(const color &mat, const color &light, double shininess, 
 *							const dvec3 &r, const dvec3 &v)
 * @brief	Computes specular color produce by a single light at a single point.
 * @param	mat		 	Material.
 * @param	lightSpecular	 	The light's color.
 * @param	shininess	Material shininess.
 * @param	r		 	Reflection vector.
 * @param	v		 	Viewing vector.
 * @return	Specular color.
 */

color specularColor(const color &mat, const color &lightSpecular,
					double shininess,
					const dvec3 &r, const dvec3 &v) {
	/* CSE 386 - todo  */
	return mat;
}

/**
 * @fn	color totalColor(const Material &mat, const LightColor &lightColor, 
 *						const dvec3 &viewingDir, const dvec3 &normal, 
 *						const dvec3 &lightPos, const dvec3 &intersectionPt, 
 *						bool attenuationOn, const LightAttenuationParameters &ATparams)
 * @brief	Color produced by a single light at a single point.
 * @param	mat			  	Material.
 * @param	lightColor	  	The light's color.
 * @param	v	  			The v vector.
 * @param	n   		  	Normal vector.
 * @param	lightPos	  	Light position.
 * @param	intersectionPt	(x,y,z) of intersection point.
 * @param	attenuationOn 	true if attenuation is on.
 * @param	ATparams	  	Attenuation parameters.
 * @return	Color produced by a single light at a single point.
 */
 
color totalColor(const Material &mat, const LightColor &lightColor,
				const dvec3 &v, const dvec3 &n,
				const dvec3 &lightPos, const dvec3 &intersectionPt,
				bool attenuationOn, 
				const LightATParams &ATparams) {
	/* CSE 386 - todo  */
	return mat.diffuse;
}

/**
 * @fn	color PositionalLight::illuminate(const dvec3 &interceptWorldCoords, 
 *										const dvec3 &normal, const Material &material, 
 *										const Frame &eyeFrame, bool inShadow) const
 * @brief	Computes the color this light produces in RAYTRACING applications.
 * @param	interceptWorldCoords	(x, y, z) at the intercept point.
 * @param	normal				The normal vector.
 * @param	material			The object's material properties.
 * @param	eyeFrame			The coordinate frame of the camera.
 * @param	inShadow			true if the point is in a shadow.
 * @return	The color produced at the intercept point, given this light.
 */

color PositionalLight::illuminate(const dvec3& interceptWorldCoords,
									const dvec3& normal,
									const Material& material,
									const Frame& eyeFrame, bool inShadow) const {
	/* CSE 386 - todo  */
	return material.diffuse;
}

/**
 * @fn	color SpotLight::illuminate(const dvec3 &interceptWorldCoords, 
 *									const dvec3 &normal, const Material &material, 
 *									const Frame &eyeFrame, bool inShadow) const
 * @brief	Computes the color this light produces in raytracing applications.
 * @param	interceptWorldCoords				The surface properties of the intercept point.
 * @param	normal					The normal vector.
 * @param	material			The object's material properties.
 * @param	eyeFrame			The coordinate frame of the camera.
 * @param	inShadow			true if the point is in a shadow.
 * @return	The color produced at the intercept point, given this light.
 */

color SpotLight::illuminate(const dvec3 &interceptWorldCoords,
							const dvec3 &normal,
							const Material &material,
							const Frame &eyeFrame, bool inShadow) const {
	/* CSE 386 - todo  */
	return material.diffuse;
}

/**
* @fn	void setDir (double dx, double dy, double dz)
* @brief	Sets the direction of the spotlight.
* @param	dx		x component of the direction
* @param	dy		y component of the direction
* @param	dz		z component of the direction
*/

void SpotLight::setDir(double dx, double dy, double dz) {
	spotDir = glm::normalize(dvec3(dx, dy, dz));
}

/**
* @fn	bool inCone(const dvec3& spotPos, const dvec3& spotDir, double spotFOV, const dvec3& intercept)
* @brief	Determines if an intercept point falls within a spotlight's cone.
* @param	spotPos		where the spotlight is positioned
* @param	spotDir		normalized direction of spotlight's pointing direction
* @param	spotFOV		spotlight's field of view, which is 2X of the angle from the viewing axis
* @param	intercept	the position of the intercept.
*/

bool inCone(const dvec3& spotPos, const dvec3& spotDir, double spotFOV, const dvec3& intercept) {
	/* CSE 386 - todo  */
	return false;
}

/**
* @fn	bool inShadow(const dvec3& lightPos, const dvec3& intercept, const dvec3& normal, const vector<VisibleIShapePtr> objects)
* @brief	Determines if an intercept point falls in a shadow.
* @param	lightPos		where the spotlight is positioned
* @param	intercept	the position of the intercept.
* @param	normal		the normal vector at the intercept point
* @param	objects		the collection of opaque objects in the scene
*/

bool inShadow(const dvec3& lightPos, const dvec3& intercept, const dvec3& normal, const vector<VisibleIShapePtr>& objects) {
	/* CSE 386 - todo  */
	return false;
}
