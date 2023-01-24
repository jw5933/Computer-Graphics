
rooms.image2D = function() {

lib2D();

description = `An interactable replica<br>of a Murakami flower...
<br>Feel free to change the colours in 'edit render' to make your own :D`;

code = {
'explanation': `
  S.html(\`
     A 2D canvas lets you create paths.
     <p>
     You can either
     draw <i>strokes</i> along those paths or else
     create solid shapes by <i>filling</i> those paths.
  \`);
`,
init: `
  S.x = 400;
  S.y = 400;
`,
assets: `
  S.line = (ax,ay,bx,by) => {
     S.context.beginPath();
     S.context.moveTo(ax,ay);
     S.context.lineTo(bx,by);
     S.context.stroke();
  }

  S.rect = (x,y,w,h) => {
     S.context.beginPath();
     S.context.rect(x,y,w,h);

     S.context.strokeStyle = 'white';
     S.context.stroke();

     if (S.isSpace) {
        S.context.fillStyle = 'gray';
        S.context.fill();
     }
  }


  //testing
  S.flower = (x, y, radius, petals, colours) => {
    for (var n = 1; n <= petals; n++) {
      S.context.fillStyle = colours[n%colours.length];
      S.context.beginPath();
      S.context.moveTo(x, y);
      let t1 = ((Math.PI * 2) / petals) * (n + 1);
      let t2 = ((Math.PI * 2) / petals) * (n);
      
      let x1 = (radius * Math.sin(t1)) + x;
      let y1 = (radius * Math.cos(t1)) + y;
      let x2 = (radius * Math.sin(t2)) + x;
      let y2 = (radius * Math.cos(t2)) + y;
      
      S.context.bezierCurveTo(x1, y1, x2, y2, 400+t1, 400+t2);
      
      S.context.stroke();
      S.context.fill();
      S.context.closePath();
    }

    //center
    // S.context.fillStyle = 'white';
    // S.context.beginPath();
    // S.context.arc(x, y, radius/5, 0, 2 * Math.PI);
    // S.context.stroke();
    // S.context.fill();
  }

  S.muraFlower = (x, y, radius, petals, colours = ['white', 'pink'], faceColour = 'yellow', mouthColour = 'red') => {
    // edit colours in 'render' to create your own flower!
    
    let c = S.context;
    c.strokeStyle = 'black';
    c.lineWidth = 5;

    for (var n = 1; n <= petals; n++) {
      c.fillStyle = colours[n%colours.length];
      c.beginPath();
      
      let centeringOffset = (Math.PI * 2 / petals)/2;
      let t1 = ((Math.PI * 2) / petals) * (n) - centeringOffset;
      let t2 = ((Math.PI * 2) / petals) * (n+1) - centeringOffset;

      let petalRadius = radius*4/6;
      //start
      let ax = petalRadius * Math.sin(t1) + x;
      let ay = petalRadius * Math.cos(t1) + y;

      //end
      let bx = petalRadius * Math.sin(t2) + x;
      let by = petalRadius * Math.cos(t2) + y;
      
      //petal curve
      let x1 = (radius * Math.sin(t1)) + x;
      let y1 = (radius * Math.cos(t1)) + y;
      let x2 = (radius * Math.sin(t2)) + x;
      let y2 = (radius * Math.cos(t2)) + y;

      c.moveTo(x, y);
      c.lineTo(ax,ay);
      c.bezierCurveTo(x1, y1, x2, y2, bx, by);
      c.moveTo(x, y);
      c.lineTo(bx,by);

      c.fill();
      c.stroke();
      
      c.closePath();
    }

    //center (face)
    const centerSize = radius*3/7;
    c.fillStyle = faceColour;
    c.beginPath();
    c.arc(x, y, centerSize, 0, 2 * Math.PI);
    c.fill();
    c.stroke();
    c.closePath();


    //eyes
    c.fillStyle = 'black';
    
    let eyeLevel = y-(centerSize)*(3/6);
    let height = radius*(1/14);
    let width = radius*(1/20);
    const distanceToEye = centerSize * 2/5;

    c.beginPath();
    c.ellipse(x - distanceToEye, eyeLevel, height, width, Math.PI/1.5, 0, 2 * Math.PI);
    c.ellipse(x + distanceToEye, eyeLevel, height, width, Math.PI/2.5, 0, 2 * Math.PI);
    c.fill();
    c.closePath();


    //eye highlights
    c.fillStyle = 'white';
    let offset = height/2.2;
    let heightDivisor = 2.7;
    let widthDivisor = 3;

    c.beginPath();
    c.ellipse(x - distanceToEye, eyeLevel - offset, height/heightDivisor, width/widthDivisor, Math.PI/1.5, 0, 2 * Math.PI);
    c.ellipse(x + distanceToEye - width/3, eyeLevel - offset, height/heightDivisor, width/widthDivisor, Math.PI/2.5, 0, 2 * Math.PI);
    c.fill();
    c.closePath();

    heightDivisor = 3.2;
    widthDivisor = 3.5;
    c.beginPath();
    c.ellipse(x - distanceToEye, eyeLevel + offset, height/heightDivisor, width/widthDivisor, Math.PI/1.5, 0, 2 * Math.PI);
    c.ellipse(x + distanceToEye + width/2.5, eyeLevel + offset, height/heightDivisor, width/widthDivisor, Math.PI/2.5, 0, 2 * Math.PI);
    c.fill();
    c.closePath();


    //happy mouth
    c.beginPath();
    const size = (centerSize) * 5/7;
    const mouthLevel = y - centerSize*(1/12);
    let w = size;
    let h = centerSize*(1/6);
    c.ellipse(x, mouthLevel, w, h, Math.PI, 0, Math.PI);
    c.fill();
    
    h = centerSize*(6/7);
    c.ellipse(x, mouthLevel, w, h, 0, 0, Math.PI);
    c.fillStyle = mouthColour;
    c.fill();
    c.stroke();
    c.closePath();
  }
`,
render: `
// edit colours to create your own flower!
  let bgColour = 'rgb(230, 245, 240)';
  let petalColours = ['white', 'rgb(242, 150, 242)', 'rgb(150, 194, 242)'];
  let faceColour = 'rgb(249, 255, 161)';
  let mouthColour = 'rgb(247, 103, 87)';
  
  
  let t = Math.sin(time);
  let c = S.context;
  
  c.fillStyle = bgColour;
  c.fillRect(0,0, 800, 800);

  const numOfPetals = 12;
  S.muraFlower(S.x, S.y, 300 /*+ 50*t*/, numOfPetals, petalColours, faceColour, mouthColour);
`,
events: `
  onDrag = (x,y) => {
     S.x = x;
     S.y = y;
  }
  onKeyPress   = key => S.isSpace = key == 32;
  onKeyRelease = key => S.isSpace = false;
`
};

}

