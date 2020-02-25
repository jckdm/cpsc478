///// OCTAHEDRON DEFINTION
/////
///// Octahedron is defined to be centered at the origin of the coordinate reference system.
///// Octahedron size is assumed to be 2.0 x 2.0 x 2.0 .
function Octahedron () {

	this.name = "octahedron";

	// vertices definition
	this.vertices = new Float32Array([
		-1.0,  0.0,  1.0,
		 1.0,  0.0,  1.0,
		 0.0,  1.0,  0.0,
		 1.0,  0.0, -1.0,
		-1.0,  0.0, -1.0,
		 0.0, -1.0,  0.0
	]);

	// triangles definition
	this.triangleIndices = new Uint16Array([
		0, 1, 2, // front
		0, 4, 2, // left
		1, 3, 2, // right
		4, 3, 2, // back
		0, 1, 5, // front-bottom
		0, 4, 5, // left-bottom
		1, 3, 5, // right-bottom
		4, 3, 5, // back-bottom
	]);

	this.numVertices = this.vertices.length/3;
	this.numTriangles = this.triangleIndices.length/3;

}
