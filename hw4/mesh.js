
rooms.mesh = function() {

lib3D2();

description = `<b>Scene 2</b>
            3D mesh made of triangle strip
            <br>
            and triangles.
	       <br><input type=range id=rate> rate
          <br><input type=range id=resolution min = 5 value = 20> resolution
	       `;

code = {
'init':`
























   // A SQUARE IS A TRIANGLE MESH WITH JUST TWO TRIANGLES

   S.squareMesh = [ -1, 1, 0,  0,0,1,  0,1,
                     1, 1, 0,  0,0,1,  1,1,
		    -1,-1, 0,  0,0,1,  0,0,
		     1,-1, 0,  0,0,1,  1,0 ];

   // GLUE TOGETHER TWO MESHES TO CREATE A SINGLE MESH

   let glueMeshes = (a,b) => {
      let mesh = a.slice();
      mesh.push(a.slice(a.length - S.VERTEX_SIZE, a.length));
      mesh.push(b.slice(0, S.VERTEX_SIZE));
      mesh.push(b);
      return mesh.flat();
   }

   // GIVEN A FUNCTION THAT MAPS (u,v) TO point AND normal,
   // AND GIVEN A MESH RESOLUTION, CREATE A PARAMETRIC MESH

   let uvMesh = (f,nu,nv) => {
      let mesh = [];
      for (let iv = 0 ; iv < nv ; iv++) {
         let v = iv / nv;
	      let strip = [];
         for (let iu = 0 ; iu <= nu ; iu++) {
	         let u = iu / nu;
            strip = strip.concat(f(u,v));
            strip = strip.concat(f(u,v+1/nv));
	      }
	      mesh = glueMeshes(mesh, strip);
      }
      return mesh;
   }

   // CREATE A UNIT SPHERE PARAMETRIC MESH
   S.sphereMesh = uvMesh((u,v) => {
      let theta = 2 * Math.PI * u;
      let phi = Math.PI * v - Math.PI/2;
      let cu = Math.cos(theta);
      let su = Math.sin(theta);
      let cv = Math.cos(phi);
      let sv = Math.sin(phi);
      return [cu * cv, su * cv, sv,
              cu * cv, su * cv, sv,
	      u, v];
   }, 20, 10);

   // CREATE A UNIT TUBE PARAMETRIC MESH
   S.tubeMesh = (res) => uvMesh((u,v) => {
      let theta = 2 * Math.PI * u;
      let phi = 2 * Math.PI * v;
      let cu = Math.cos(theta);
      let su = Math.sin(theta);
      return [cu, su, 2 * v -1,
              cu, su, 0, u, v];
   }, res, 2);

   // CREATE A UNIT CONE PARAMETRIC MESH
   S.openConeMesh = (res) => 
      uvMesh((u,v) => {
      let theta = 2 * Math.PI * u;
      let phi = 2 * Math.PI * v;
      let cu = Math.cos(theta);
      let su = Math.sin(theta);
      return [v * cu, v * su, 2*v-1,
              cu, su, 0, u, v];
   }, res, 2);

   // CREATE A UNIT DISK PARAMETRIC MESH
   S.diskMesh = (res) => uvMesh((u,v) => {
      let theta = 2 * Math.PI * u;
      let phi = 2 * Math.PI * v;
      let cu = Math.cos(theta);
      let su = Math.sin(theta);
      return [v * cu, v * su, 0,
              0, 0, 1, u, v];
   }, res, 2);


   // TRANSFORM A MESH BY A MATRIX ON THE CPU

   let transformMesh = (mesh, matrix) => {
      let result = [];
      let IMT = matrixTranspose(matrixInverse(matrix));
      for (let n = 0 ; n < mesh.length ; n += S.VERTEX_SIZE) {
         let V = mesh.slice(n, n + S.VERTEX_SIZE);
	 let P  = V.slice(0, 3);
	 let N  = V.slice(3, 6);
	 let UV = V.slice(6, 8);
	 P = matrixTransform(matrix, [P[0], P[1], P[2], 1]);
	 N = matrixTransform(IMT,    [N[0], N[1], N[2], 0]);
         result.push(P[0],P[1],P[2], N[0],N[1],N[2], UV);
      }
      return result.flat();
   }

   // A CUBE MESH IS SIX TRANSFORMED SQUARE MESHES GLUED TOGETHER

   let face0 = transformMesh(S.squareMesh, matrixTranslate([0,0,1]));
   let face1 = transformMesh(face0,        matrixRotx( Math.PI/2));
   let face2 = transformMesh(face0,        matrixRotx( Math.PI  ));
   let face3 = transformMesh(face0,        matrixRotx(-Math.PI/2));
   let face4 = transformMesh(face0,        matrixRoty(-Math.PI/2));
   let face5 = transformMesh(face0,        matrixRoty( Math.PI/2));
   S.cubeMesh = glueMeshes(face0,
                glueMeshes(face1,
                glueMeshes(face2,
                glueMeshes(face3,
                glueMeshes(face4,
		           face5)))));
   
   //disk ends
   let end = (res) => transformMesh(S.diskMesh(res), matrixTranslate([0,0,1]));
   let endOpposite = (res) => transformMesh(end(res), matrixRotx(Math.PI));

   //cone mesh
   S.coneMesh = (res) => glueMeshes(S.openConeMesh(res), end(res));

   //cylinder mesh
   S.cylinderMesh = (res) => glueMeshes(S.tubeMesh(res), glueMeshes(end(res), endOpposite(res)));

   // DRAW A SINGLE MESH. WE STILL NEED TO ADD MATERIAL PROPERTIES!
   S.drawMesh = (mesh, matrix) => {
      let gl = S.gl;
      S.setUniform('Matrix4fv', 'uMatrix', false, matrix);
      S.setUniform('Matrix4fv', 'uInvMatrix', false, matrixInverse(matrix));
      S.gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(mesh), gl.STATIC_DRAW);
      S.gl.drawArrays(S.gl.TRIANGLE_STRIP, 0, mesh.length / S.VERTEX_SIZE);
   }

   S.triangleSquare = [
      0, 0, 0,    0, -1, 0,    0, 0,
      1, 0, 0,    0, -1, 0,    1, 0,
      .5, 0, .5,    0, -1, 0,    0, 1,

      0, 0, 0,    0, -1, 0,    0, 0,
      0, 0, 1,    0, -1, 0,    1, 0,
      .5, 0, .5,    0, -1, 0,    0, 1,

      1, 0, 1,    0, -1, 0,    0, 0,
      1, 0, 0,    0, -1, 0,    1, 0,
      .5, 0, .5,    0, -1, 0,    0, 1,

      1, 0, 1,    0, -1, 0,    0, 0,
      0, 0, 1,    0, -1, 0,    1, 0,
      .5, 0, .5,    0, -1, 0,    0, 1,
   ];

   S.pyramid = [
      0, 0, 0,    -0.5, 0, -0.5,    0, 0,
      1, 0, 0,    0.5, 0, -0.5,    1, 0,
      .5, 1, .5,    0, 0, -1,    0, 1,

      0, 0, 0,    -0.5, 0, -0.5,    0, 0,
      0, 0, 1,    -0.5, 0, 0.5,    1, 0,
      .5, 1, .5,    -1, 0, 0,    0, 1,

      1, 0, 1,    0.5, 0, 0.5,    0, 0,
      0, 0, 1,    -0.5, 0, 0.5,    1, 0,
      .5, 1, .5,    0, 0, 1,    0, 1,

      1, 0, 1,    0.5, 0, 0.5,    0, 0,
      1, 0, 0,    0.5, 0, -0.5,    1, 0,
      .5, 1, .5,    1, 0, 0,    0, 1,

      0, 0, 0,    0, -1, 0,    0, 0,
      1, 0, 0,    0, -1, 0,    1, 0,
      .5, 0, .5,    0, -1, 0,    0, 1,

      0, 0, 0,    0, -1, 0,    0, 0,
      0, 0, 1,    0, -1, 0,    1, 0,
      .5, 0, .5,    0, -1, 0,    0, 1,

      1, 0, 1,    0, -1, 0,    0, 0,
      1, 0, 0,    0, -1, 0,    1, 0,
      .5, 0, .5,    0, -1, 0,    0, 1,

      1, 0, 1,    0, -1, 0,    0, 0,
      0, 0, 1,    0, -1, 0,    1, 0,
      .5, 0, .5,    0, -1, 0,    0, 1,
   ];
   S.octahedron = [
      0, 0, 0,    -0.5, 0, -0.5,    0, 0,
      1, 0, 0,    0.5, 0, -0.5,    1, 0,
      .5, 1, .5,    0, 0, -1,    0, 1,

      0, 0, 0,    -0.5, 0, -0.5,    0, 0,
      0, 0, 1,    -0.5, 0, 0.5,    1, 0,
      .5, 1, .5,    -1, 0, 0,    0, 1,

      1, 0, 1,    0.5, 0, 0.5,    0, 0,
      0, 0, 1,    -0.5, 0, 0.5,    1, 0,
      .5, 1, .5,    0, 0, 1,    0, 1,

      1, 0, 1,    0.5, 0, 0.5,    0, 0,
      1, 0, 0,    0.5, 0, -0.5,    1, 0,
      .5, 1, .5,    1, 0, 0,    0, 1,
      
      0, 0, 0,    -0.5, 0, -0.5,    0, 0,
      1, 0, 0,    0.5, 0, -0.5,    1, 0,
      .5, -1, .5,    0, 0, -1,    0, 1,

      0, 0, 0,    -0.5, 0, -0.5,    0, 0,
      0, 0, 1,    -0.5, 0, 0.5,    1, 0,
      .5, -1, .5,    -1, 0, 0,    0, 1,

      1, 0, 1,    0.5, 0, 0.5,    0, 0,
      0, 0, 1,    -0.5, 0, 0.5,    1, 0,
      .5, -1, .5,    0, 0, 1,    0, 1,

      1, 0, 1,    0.5, 0, 0.5,    0, 0,
      1, 0, 0,    0.5, 0, -0.5,    1, 0,
      .5, -1, .5,    1, 0, 0,    0, 1,
   ];

`,
fragment: `
S.setFragmentShader(\`
   varying vec3 vPos, vNor;
   void main() {
      float c = .2 + .8 * max(0.,dot(vNor,vec3(.57)));
      gl_FragColor = vec4(c,c,c,1.);
   }
\`);
`,
vertex: `
S.setVertexShader(\`

   attribute vec3 aPos, aNor;
   varying   vec3 vPos, vNor;
   uniform   mat4 uMatrix, uInvMatrix, uProject;

   void main() {
      vec4 pos = uProject * uMatrix * vec4(aPos, 1.);
      vec4 nor = vec4(aNor, 0.) * uInvMatrix;
      vPos = pos.xyz;
      vNor = normalize(nor.xyz);
      gl_Position = pos * vec4(1.,1.,-.01,1.);
   }
\`)
`,
render: `

   // SET THE PROJECTION MATRIX BASED ON CAMERA FOCAL LENGTH

   let fl = 5.0;
   S.setUniform('Matrix4fv', 'uProject', false,
      [1,0,0,0, 0,1,0,0, 0,0,1,-1/fl, 0,0,0,1]);

   let m = new Matrix();

   // GET VALUES FROM THE HTML SLIDERS

   let T = 2 * time * rate.value / 100;
   let res = resolution.value;

   // RENDER THE SCENE

   m.identity();
   m.save();

   // HEAD
	   m.save();
         m.scale(.2,.2,.2);
         m.roty(.5 * T);
         S.drawMesh(S.coneMesh(res), m.get());
      m.restore();
      m.save();
         m.scale(.2,.2,.2);
         m.translate(1, 1, 1);
         m.roty(.5 * T);
         m.rotz(.2 * T);
         S.drawMesh(S.octahedron, m.get());
      m.restore();
      m.restore();

`,
events: `
   ;
`
};

}

