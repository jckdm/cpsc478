// adapted from http://stackoverflow.com/questions/3914203/javascript-filter-image-color
function createCanvas(image)
{
  // create a new canvas element
  var myCanvas = document.createElement("canvas");
  var myCanvasContext = myCanvas.getContext("2d");

  var imgWidth=image.width;
  var imgHeight=image.height;

  // set the width and height to show two copies of the image
  myCanvas.width= 2*imgWidth+ 10;
  myCanvas.height = imgHeight;

  // draw the image
  myCanvasContext.drawImage(image,0,0);

  // get all the input and output image data into arrays
  var imageData = myCanvasContext.getImageData(0,0, imgWidth, imgHeight);
  var imoutData = myCanvasContext.getImageData(0,0, imgWidth, imgHeight);

  // go through it all...
  for (j=0; j<imageData.width; j++) {
    for (i=0; i<imageData.height; i++) {
       var red = 0;
       var green = 0;
       var blue = 0;
       var num = 0;

      // ADAPTED (and translated into JS) FROM STAFF SOLUTION OF CS50 FALL 19 PSET "filter/blur"
      for (var r = Math.max(0, i - 2); r <= Math.min(imageData.height - 1, i + 2); r++) {
          for (var c = Math.max(0, j - 2); c <= Math.min(imageData.width - 1, j + 2); c++) {

            var index = (r*4)*imageData.width+(c*4);
            red += imageData.data[index];
            green += imageData.data[index+1];
            blue += imageData.data[index+2];
            num++;
          }
       }
       imoutData.data[index]   = Math.round(red / num);
       imoutData.data[index+1] = Math.round(green / num);
       imoutData.data[index+2] = Math.round(blue / num);
     }
   }

   // put the image data back into the canvas
   myCanvasContext.putImageData(imoutData, imageData.width+10,0,0,0, imageData.width, imageData.height);

  // append it to the body
  document.body.appendChild(myCanvas);
}

function loadImage() {
  var img = new Image();
  img.onload = function () {
    createCanvas(img);
  }
  img.src = document.getElementById("imagefilename").files[0].name;
}
