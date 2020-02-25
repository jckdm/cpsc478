function doStuff() {
  var c = document.getElementById("myCanvas");
  var context = c.getContext("2d");
  var imgWidth = 500;
  var imgHeight = 500;
  var probability = document.getElementById("probability").value;

  var choices = prob(parseInt(probability * 10));

  var onecolor = [255, 255, 255];
  onecolor[0] = document.getElementById("oneRed").value;
  onecolor[1] = document.getElementById("oneGreen").value;
  onecolor[2] = document.getElementById("oneBlue").value;

  var twocolor = [255, 255, 255];
  twocolor[0] = document.getElementById("twoRed").value;
  twocolor[1] = document.getElementById("twoGreen").value;
  twocolor[2] = document.getElementById("twoBlue").value;

  var imgData = context.createImageData(imgWidth, imgHeight);

  for (jj = 0; jj < imgWidth; jj++) {
    var qq = jj*4;
    for(ii = 0; ii < imgHeight; ii += 1) {
      var pp = (ii*imgWidth*4)+qq;
      var choice = choices[Math.floor(Math.random() * 10)]

      if (choice == 1) {
        imgData.data[pp]   = onecolor[0];
        imgData.data[pp+1] = onecolor[1];
        imgData.data[pp+2] = onecolor[2];
        imgData.data[pp+3] = 255;
      }
      else {
        imgData.data[pp]   = twocolor[0];
        imgData.data[pp+1] = twocolor[1];
        imgData.data[pp+2] = twocolor[2];
        imgData.data[pp+3] = 255;
      }
    }
  }
  context.putImageData(imgData, 20, 20);
}

function prob(p) {
  var arr = [];
  arr.length = 10;

  for (var i = 0; i < 10; i++) {
    if (i < p) { arr[i] = 1; }
    else { arr[i] = 2; }
  }
  return arr;
}
