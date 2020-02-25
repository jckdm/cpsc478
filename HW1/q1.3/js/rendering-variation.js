var shaderProgram  = null;
var vertexBuffer = null;
var vertexColorBuffer = null;
var aPositionIndex = -1;
var aVertexColor = -1;

///// Initialize the data buffer to pass to the rendering pipeline
///// the geometry and its attributes.
function initBuffers(gl)
{
  var triangles = document.getElementById('tri').value;
  var l = triangles * 9;

  triangleVertices = new Float32Array(l+9);

  var start = [0, 0, 0];
  // plus one for final triangle
  for (var i = 0; i < triangles+1; i++) {

    var j = i * 9;

    triangleVertices[j    ] = start[0];
    triangleVertices[j + 1] = start[1];
    triangleVertices[j + 2] = start[2];

    start[0] = Math.cos(i * (Math.PI * 2/triangles));
    start[1] = Math.sin(i * (Math.PI * 2/triangles));
    start[2] = 0;

    triangleVertices[j + 3] = start[0];
    triangleVertices[j + 4] = start[1];
    triangleVertices[j + 5] = start[2];

    triangleVertices[j + 6] = 0;
    triangleVertices[j + 7] = 0;
    triangleVertices[j + 8] = 0;
  }

  var onecolor = [1.0, 1.0, 1.0];
  onecolor[0] = document.getElementById("oneRed").value;
  onecolor[1] = document.getElementById("oneGreen").value;
  onecolor[2] = document.getElementById("oneBlue").value;

  var twocolor = [1.0, 1.0, 1.0];
  twocolor[0] = document.getElementById("twoRed").value;
  twocolor[1] = document.getElementById("twoGreen").value;
  twocolor[2] = document.getElementById("twoBlue").value;

  var triangleVerticesColor = new Float32Array(l+9), color;

  if (document.getElementById("circus").checked) {
    for (var i = 0; i < triangles * 3 + 3; i++) {
      if (i % 2 != 0) { color = onecolor; }
      else { color = twocolor; }
      for (var j = 0; j < 3; j++) {
        triangleVerticesColor[i * 9 + j * 3    ] = color[0];
        triangleVerticesColor[i * 9 + j * 3 + 1] = color[1];
        triangleVerticesColor[i * 9 + j * 3 + 2] = color[2];
      }
    }
  }
  // RAMP
  else {
    for (var i = 0; i < triangles * 3 + 3; i++) {
      for (var j = 0; j < 3; j++) {
        if (j != 2) { color = onecolor; }
        else { color = twocolor; }
        triangleVerticesColor[i * 9 + j * 3    ] = color[0];
        triangleVerticesColor[i * 9 + j * 3 + 1] = color[1];
        triangleVerticesColor[i * 9 + j * 3 + 2] = color[2];
      }
    }
  }

  vertexBuffer = gl.createBuffer();
  gl.bindBuffer(gl.ARRAY_BUFFER, vertexBuffer);
  gl.bufferData(gl.ARRAY_BUFFER, triangleVertices, gl.STATIC_DRAW);
  gl.bindBuffer(gl.ARRAY_BUFFER, null);

  vertexColorBuffer = gl.createBuffer();
  gl.bindBuffer(gl.ARRAY_BUFFER, vertexColorBuffer);
  gl.bufferData(gl.ARRAY_BUFFER, triangleVerticesColor, gl.STATIC_DRAW);
  gl.bindBuffer(gl.ARRAY_BUFFER, null);
}

///// Define and compile a very simple shader.
function initShaders(gl)
{
  var vertexShaderSource = "\
  attribute vec3 a_position;                  \n\
  attribute vec3 a_color;                     \n\
  varying vec3 vertexcolor;                   \n\
  void main(void)                             \n\
  {                                           \n\
    vertexcolor = a_color;                  \n\
    gl_Position = vec4(a_position, 1.0);    \n\
  }                                           \n\
  ";

  var fragmentShaderSource = "\
  precision highp float;                      \n\
  varying vec3 vertexcolor;                   \n\
  void main(void)                             \n\
  {                                           \n\
    gl_FragColor = vec4(vertexcolor, 1.0);  \n\
  }                                           \n\
  ";

  // create the vertex shader
  var vertexShader = gl.createShader(gl.VERTEX_SHADER);
  gl.shaderSource(vertexShader, vertexShaderSource);
  gl.compileShader(vertexShader);

  // create the fragment shader
  var fragmentShader = gl.createShader(gl.FRAGMENT_SHADER);
  gl.shaderSource(fragmentShader, fragmentShaderSource);
  gl.compileShader(fragmentShader);

  // Create the shader program
  shaderProgram = gl.createProgram();
  gl.attachShader(shaderProgram, vertexShader);
  gl.attachShader(shaderProgram, fragmentShader);
  gl.linkProgram(shaderProgram);

  // If creating the shader program failed, we show compilation and linking errors.
  if (!gl.getProgramParameter(shaderProgram, gl.LINK_STATUS))
  {
    alert("Unable to initialize the shader program.");
    var str = "";
    str += "VS:\n" + gl.getShaderInfoLog(vertexShader) + "\n\n";
    str += "FS:\n" + gl.getShaderInfoLog(fragmentShader) + "\n\n";
    str += "PROG:\n" + gl.getProgramInfoLog(shaderProgram);
    alert(str);
  }
}

///// Draw the given triangle interpolating vertices color.
function renderTriangle(gl)
{
  // Clear the framebuffer of the rendering context
  gl.clearColor(0.0, 0.0, 0.0, 1.0);
  gl.clear(gl.COLOR_BUFFER_BIT);

  // enable the current shader program
  gl.useProgram(shaderProgram);

  // connect the buffer containing the vertices of the triangle with the position attribute
  gl.bindBuffer(gl.ARRAY_BUFFER, vertexBuffer);
  aPositionIndex = gl.getAttribLocation(shaderProgram, "a_position");
  gl.enableVertexAttribArray(aPositionIndex);
  gl.vertexAttribPointer(aPositionIndex, 3, gl.FLOAT, false, 0, 0);

  // connect the buffer containing the color of each vertex with the color attribute
  gl.bindBuffer(gl.ARRAY_BUFFER, vertexColorBuffer);
  aVertexColor = gl.getAttribLocation(shaderProgram, "a_color");
  gl.enableVertexAttribArray(aVertexColor);
  gl.vertexAttribPointer(aVertexColor, 3, gl.FLOAT, false, 0, 0);

  // start to draw (!)
  var triangles = document.getElementById('tri').value;
  gl.drawArrays(gl.TRIANGLES, 0, (triangles*3)+3);

  // disable the current shading program
  gl.useProgram(null);
}
