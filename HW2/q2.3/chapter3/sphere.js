function midpoint(arr1, arr2) {
	var arr = new Float32Array(3);

	arr[0] = (arr1[0] + arr2[0]) / 2;
	arr[1] = (arr1[1] + arr2[1]) / 2;
	arr[2] = (arr1[2] + arr2[2]) / 2;

	return arr;
}

function Sphere(n) {
	this.name = "sphere";

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

	for (var i = 0; i < n; i++) {
		for (var j = 0; j < this.numTriangles; j++) {

			var ind1 = 3 * this.triangleIndices[3*j    ];
			var ind2 = 3 * this.triangleIndices[3*j + 1];
			var ind3 = 3 * this.triangleIndices[3*j + 2];

			var vert1 = this.vertices.slice(ind1, ind1 + 3);
			var vert2 = this.vertices.slice(ind2, ind2 + 3);
			var vert3 = this.vertices.slice(ind3, ind3 + 3);

			var mid1 = midpoint(vert1, vert2);
			var mid2 = midpoint(vert1, vert3);
			var mid3 = midpoint(vert2, vert3);

			var tempV = Array.from(this.vertices);
			var tempI = Array.from(this.triangleIndices);

			tempV.push(mid1[0], mid1[1], mid1[2], mid2[0], mid2[1], mid2[2], mid3[0], mid3[1], mid3[2]);

			var l = tempV.length / 3;

			tempI.push(ind1/3, l-2, l-3,ind2/3, l-1, l-3, ind3/3, l-2, l-1, l-1, l-2, l-3);

			var temp32 = new Float32Array(tempV);
			var temp16 = new Uint16Array(tempI);

			this.vertices = temp32;
			this.triangleIndices = temp16;
		}
		this.numVertices = this.vertices.length/3;
		this.numTriangles = this.triangleIndices.length/3;
	}

	for (var i = 0; i < this.numVertices; i++) {
		var x = this.vertices[3*i   ];
		var y = this.vertices[3*i + 1];
		var z = this.vertices[3*i + 2];

		var scalar = Math.sqrt(1 / (Math.pow(x, 2) + Math.pow(y, 2) + Math.pow(z, 2)));

		this.vertices[3*i    ] = x * scalar;
		this.vertices[3*i + 1] = y * scalar;
		this.vertices[3*i + 2] = z * scalar;
	}
}
