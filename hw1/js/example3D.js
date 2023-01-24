
rooms.image3D = function() {

lib3D();

description = 'An attempt at creating a shoreline.<br>It gets crazy as time goes on.';

code = {
'explanation': `
   S.html(\`
      Most of the work happens in a fragment shader.
      <p>
      Input to the fragment shader is x,y and time: <code>uPos, uTime</code>
      <p>
      We can also interact by adding information about the cursor: <code>uX,uY</code>
      <p>
      Output at each fragment is: red,green,blue,alpha
   \`);
`,

vertex: `
S.setVertexShader(\`

   attribute vec3 aPos;
   varying   vec3 vPos;

   void main() {
      vPos = aPos;
      gl_Position = vec4(aPos, 1.);
   }

\`)
`,

fragment: `

S.setFragmentShader(\`

   uniform float uTime, uSpace, uX, uY;
   varying vec3 vPos;

   float turbulence(vec3 p) {
      float t = 0., f = 1.;
      for (int i = 0 ; i < 5 ; i++) {
         t += abs(noise(f * p)) / f;
         f *= 3.;
      }
      return t;
   }

   float fractal(vec3 p) {
      float t = 0., f = 1.;
      for (int i = 0 ; i < 10 ; i++) {
         t += noise(f * p) / f;
         f *= 2.;
      }
      return t;
   }
   
   
   vec3 shoreColour(float y) {
      vec3 water = vec3(.13,.85,.67);
      water = mix(water, vec3(0.89,0.86,0.60), max(y+.8,1.));
      vec3 shore = mix(water, vec3(.2,.5,.9), max(2.2*y, 0.));
      return shore;
   }

   void main() {
      vec3 p0 = vPos + vec3(0.2, 0.,.03*uTime);
      vec3 p1 = vPos + vec3(0.,.5,.03*uTime);
      vec3 color = shoreColour(p0.y + turbulence(mix(p0,p1, cos(.5*uTime)) * fractal(mix(p0,p1, cos(.7*uTime))) ));

      gl_FragColor = vec4(sqrt(color), 1.);
   }

\`)

`,

render: `
   S.setUniform('1f', 'uTime', time);
   S.gl.drawArrays(S.gl.TRIANGLE_STRIP, 0, 4);
`,

events: `
   onDrag = (x,y) => {
      S.setUniform('1f', 'uX', x);
      S.setUniform('1f', 'uY', y);
   }
   onKeyPress  =k=>S.setUniform('1f','uSpace',k==32);
   onKeyRelease=k=>S.setUniform('1f','uSpace',false);
`
}

}

