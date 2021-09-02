/****************************************************
 * 2016-2021 Eric Bachmann and Mike Zmuda
 * All Rights Reserved.
 * NOTICE:
 * Dissemination of this information or reproduction
 * of this material is prohibited unless prior written
 * permission is granted.
 ****************************************************/

#include <ctime>
#include <vector>
#include "defs.h"
#include "utilities.h"
#include "framebuffer.h"
#include "colorandmaterials.h"
#include "rasterization.h"
#include "io.h"

FrameBuffer frameBuffer(WINDOW_WIDTH, WINDOW_HEIGHT);

const int SZ = 51;
const int SZ2 = SZ / 2;

void closedSquare(int x, int y, color C) {
}

void closedSquare(const ivec2 &centerPt, color C) {
}

void openSquare(const ivec2 &centerPt, color C) {
}

void render() {
	frameBuffer.clearColorAndDepthBuffers();

	drawLine(frameBuffer, 0, 0, 100, 100, red);
	drawLine(frameBuffer, 100, 100, 200, 100, green);

	closedSquare(100, 150, red);
	closedSquare(ivec2(200, 150), green);
	openSquare(ivec2(300, 150), blue);

	frameBuffer.showColorBuffer();
}

void resize(int width, int height) {
	frameBuffer.setFrameBufferSize(width, height);
	glutPostRedisplay();
}
int main(int argc, char *argv[]) {
 //   graphicsInit(argc, argv, __FILE__);
 //       
	//glutDisplayFunc(render);
	//glutReshapeFunc(resize);
	//glutKeyboardFunc(keyboardUtility);
	//glutMouseFunc(mouseUtility);

	//frameBuffer.setClearColor(white);

	//glutMainLoop();

	//double x1, y1;
	//pointOnUnitCircle(PI_2, x1, y1);
	//cout << x1 << ' ' << y1 << endl;

	//dvec2 center(3.0, 2.0);
	//dvec2 pt = pointOnCircle(center, 2.0, PI_2);
	//cout << pt << endl;

	dvec2 target(0.0, 0.0);
	cout << directionInRadians(target) << endl;

	//dvec2 target2(0.0, -2.0);
	//cout << directionInRadians(target2) << endl;

	//dvec2 targetAgain(-2.0, -2.0);
	//cout << directionInRadians(targetAgain) << endl;

	//dvec2 reference(0.0, 0.0);
	//dvec2 target3(2.0, 2.0);
	//cout << directionInRadians(reference, target3) << endl;

	//dvec2 reference2(2.0, 10.0);
	//dvec2 target4(3.0, 11.0);
	//cout << directionInRadians(reference2, target4) << endl;

	//dvec2 reference3(2.0, 2.0);
	//dvec2 target5(2.0, 0.0);
	//cout << directionInRadians(reference3, target5) << endl;

	//cout << directionInRadians(0, 0, 2, 2) << endl;
	//cout << directionInRadians(2, 10, 3, 11) << endl;
	//cout << directionInRadians(2, 2, 2, 0) << endl;

	// labtest

	return 0;
}