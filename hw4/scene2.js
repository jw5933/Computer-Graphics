
rooms.scene2 = function() {

lib3D2();

description = `<b>Scene 2</b>
               <p>
               Hierarchic 3D scene
               <br>
               with triangle meshes.
	       <p> 
          <input type=range id=pelvis_length value=25 min = 15 max = 35>  body length
          <br><input type=range id=arm_length value=30 min = 20 max = 40> arm length
	       <br><input type=range id=leg_length value=40 min = 30 max = 50> leg length
	       <br><input type=range id=rate> rate
	       `;

code = {
'init':`
























   //lighting
   S.nL = 2;
   S.jointColor = [.25,.34,.8,         0,.5,0,0,0,       2,2,2,20,      0,0,0,0];
   S.gold =    [.25,.15,.025,    0,.5,.3,.05,0,    1,.6,.1,6,     0,0,0,0];
   S.plastic = [.003,.6,.9,         0,.5,0,0,0,       2,2,2,20,      0,0,0,0];
   S.copper =  [.15,.05,.025,0,  .3,.1,.05,0,      .6,.2,.1,3,    0,0,0,0];
   S.silver =  [.1,.1,.1,0,      .1,.1,.1,0,       1,1,1,5,       0,0,0,0];

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

   S.drawMeshColor = (mesh, matrix, phong) => {
      let gl = S.gl;
      S.setUniform('Matrix4fv', 'uMatrix', false, matrix);
      S.setUniform('Matrix4fv', 'uInvMatrix', false, matrixInverse(matrix));
      S.setUniform('Matrix4fv', 'uMaterial', false, phong);
      S.gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(mesh), gl.STATIC_DRAW);
      S.gl.drawArrays(S.gl.TRIANGLE_STRIP, 0, mesh.length / S.VERTEX_SIZE);
   }

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
   const int nL = \` + S.nL + \`;
   varying vec3 vPos, vNor;
   uniform mat4 uMaterial;
   uniform vec3 uLd[nL];
   uniform vec3 uLc[nL];

   vec3 shadeSurface(vec3 P, vec3 N, mat4 M){
      vec3 ambient  = M[0].rgb;
      vec3 diffuse  = M[1].rgb;
      vec3 specular = M[2].rgb;
      float p       = M[2].a;

      // vec3 c = mix(ambient, uBgColor, .3);
      vec3 c = ambient;
      vec3 E = vec3(0.,0.,1.); //eye

      for (int l = 0 ; l < nL ; l++) {
         
         float t = -1.;

         // Lighting if not in shadow
         if (t < 0.) {
            vec3 R = 2. * dot(N, uLd[l]) * N - uLd[l];
            c += uLc[l] * (diffuse * max(0.,dot(N, uLd[l]))
                        + specular * pow(max(0., dot(R, E)), p));
         }
      }

      return c;
   }

   void main() {
      vec3 color = shadeSurface(vPos, vNor, uMaterial);
      float c = .2 + .8 * max(0.,dot(vNor,vec3(.57)));
      gl_FragColor = vec4(color,1.);
      // gl_FragColor = vec4(c,c,c,1.);
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
   
   // SEND LIGHT SOURCE DATA TO GPU
   let dot = (a,b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
   let norm = v => Math.sqrt(dot(v,v));
   let normalize = v => { let s = norm(v); return [ v[0]/s, v[1]/s, v[2]/s ]; }

   let ldData = [ normalize([1,1,1]),
                  normalize([-1,-1,-1]) ];
   S.setUniform('3fv', 'uLd', ldData.flat());
   S.setUniform('3fv', 'uLc', [ 1,1,1, .5,.3,.1 ]);

   // SET THE PROJECTION MATRIX BASED ON CAMERA FOCAL LENGTH

   let fl = 5.0;
   S.setUniform('Matrix4fv', 'uProject', false,
      [1,0,0,0, 0,1,0,0, 0,0,1,-1/fl, 0,0,0,1]);

   let m = new Matrix();

   // GET VALUES FROM THE HTML SLIDERS

   let T = 2 * time * rate.value / 100;
   let AL = .1 + .9 * arm_length.value / 100;
   let LL = .1 + .9 * leg_length.value / 100;
   let PL = .1 + .9 * pelvis_length.value /100;

   let jSize = (0.04, 0.04, 0.04);

   // RENDER THE SCENE

m.identity();
m.scale(0.9);

//pelvis
m.save();
   m.roty(Math.sin(T)*0.5);
   // m.rotz(Math.sin(T));
   // m.rotx(Math.sin(T)*0.5);
   
   //joint
   m.save();
      m.translate(-jSize, 0, 0);
      m.scale(jSize*2);
      S.drawMeshColor(S.octahedron, m.get(), S.copper);
   m.restore();

   //pelvis rigid
   m.save();
      m.scale(PL * 0.55, PL * 0.5, PL * 0.5 * 0.9);
      S.drawMeshColor(S.tubeMesh(20), m.get(), S.gold);
   m.restore();


   //legs
   for (let s = -1 ; s <= 1 ; s += 2) {
      let t = T + s * Math.PI/2;
      //hip
      m.save();
         m.translate(s* PL * 0.55, -0.1, 0);
         m.rotx(-0.5 + Math.cos(t));
         m.rotz(Math.sin(t)*0.5 + 0.3*s);
         
         //joint
         m.save();
            m.scale(jSize);
            S.drawMeshColor(S.sphereMesh, m.get(), S.jointColor);
         m.restore();

         //thigh rigid
         m.save();
            m.translate(0,-LL/2,0);
            m.scale(.05,LL/2,.05);
            S.drawMeshColor(S.sphereMesh, m.get(), S.plastic);
         m.restore();
         
         //knee
         m.save();
            m.translate(0, -LL, 0);
            m.rotx(1 + Math.sin(t));
            
            //joint
            m.save();
               m.scale(jSize);
               S.drawMeshColor(S.sphereMesh, m.get(), S.jointColor);
            m.restore();

            //shin rigid
            m.save();
               m.translate(0,-LL/2,0);
               m.scale(.05,LL/2,.05);
               S.drawMeshColor(S.sphereMesh, m.get(), S.plastic);
            m.restore();

            //ankle
            m.save();
               m.translate(0, -LL, 0);
               // m.roty(Math.sin(t));
               m.rotx(-Math.sin(t)*0.4 - Math.PI/3);
               
               //joint
               m.save();
                  m.scale(jSize);
                  S.drawMeshColor(S.sphereMesh, m.get(), S.jointColor);
               m.restore();
               
               //foot rigid
               m.save();
                  m.translate(0, -0.1, 0);
                  m.scale(0.05, 0.08, 0.03);
                  S.drawMeshColor(S.cubeMesh, m.get(), S.gold);
               m.restore();
            m.restore();
         m.restore();
      m.restore();
   }

   //waist (top half)
   m.save();
      m.translate(0, PL * 0.5, 0);
      m.rotx(Math.cos(T) * 0.5);
      m.roty(Math.cos(T));
      
      //waist joint
      m.save();
         m.scale(jSize);
         S.drawMeshColor(S.sphereMesh, m.get(), S.jointColor);
      m.restore();

      //chest rigid
      let cH = PL * 1.3 * 0.5;
      m.save();
         m.translate(0, cH, 0);
         m.scale(cH * 0.9, cH, cH * 0.6);
         S.drawMeshColor(S.sphereMesh, m.get(), S.gold);
      m.restore();

      //collarbone
      m.save();
         m.translate(0, PL * 1.3, 0);
         m.roty(Math.sin(T)*0.3);
         m.rotx(Math.cos(T)*0.7);

         //joint
         m.save();
            m.scale(jSize);
            S.drawMeshColor(S.sphereMesh, m.get(), S.jointColor);
         m.restore();

      
         //arms
         for (let s = -1 ; s <= 1 ; s += 2) {
            let t = T + s * Math.PI/2;
            
            //arm
            m.save();
               let sW = cH * 1.2;
               m.rotz(Math.sin(t)*0.3);
               m.translate(s * sW, 0, 0);
               
               //shoulder blade rigid
               m.save();
                  m.translate(-s*sW*0.5, 0, 0);
                  m.scale(sW/2, 0.015, 0.03);
                  S.drawMeshColor(S.sphereMesh, m.get(), S.gold);
               m.restore();

               //shoulder
               m.save();
                  m.rotx(Math.sin(t));
                  m.rotz(Math.sin(t)*0.3 + .6*s);
                  
                  //joint
                  m.save();
                     m.scale(jSize);
                     S.drawMeshColor(S.sphereMesh, m.get(), S.jointColor);
                  m.restore();

                  //upper arm rigid
                  m.save();
                     m.translate(0,-AL/2,0);
                     m.scale(.05, AL/2, .05);
                     S.drawMeshColor(S.sphereMesh, m.get(), S.plastic);
                  m.restore();
                  

                  //elbow
                  m.save();
                     m.translate(0,-AL,0);
                     m.rotx(-1 + .7 * Math.sin(t));
                     //joint
                     m.save();
                        m.scale(jSize);
                        S.drawMeshColor(S.sphereMesh, m.get(), S.jointColor);
                     m.restore();
                     
                     //forearm rigid
                     m.save();
                        m.translate(0,-AL/2,0);
                        m.scale(.05, AL/2, .05);
                        S.drawMeshColor(S.sphereMesh, m.get(), S.plastic);
                     m.restore();
                     
                     //wrist
                     m.save();
                        m.translate(0,-AL,0);
                        m.rotz(Math.sin(t)*0.5);
                        m.roty(Math.sin(t+0.0235));
                        
                        //joint
                        m.save();
                           m.scale(jSize);
                           S.drawMeshColor(S.sphereMesh, m.get(), S.jointColor);
                        m.restore();
                        
                        //hand rigid
                        m.save();
                           m.translate(0, -0.1, 0);
                           m.scale(0.05, 0.08, 0.03);
                           S.drawMeshColor(S.cubeMesh, m.get(), S.gold);
                        m.restore();
                     m.restore();

                  m.restore();
               m.restore();
            m.restore();
         }

         //neck rigid
         m.save();
            m.translate(0, 0.06, 0);
            m.scale(0.03, 0.06, 0.03);
            S.drawMeshColor(S.cubeMesh, m.get(), S.gold);
         m.restore();

         //head
         m.save();
            m.translate(0, 0.12, 0);
            m.roty(Math.cos(T)*0.3);
            m.rotx(Math.sin(T+0.2)*0.3 + 0.3);
            
            //joint
            m.save();
               m.scale(jSize);
               S.drawMeshColor(S.sphereMesh, m.get(), S.jointColor);
            m.restore();
            
            //head rigid
            m.save();
               m.translate(0, 0.13, 0);
               m.scale(0.1, 0.13, 0.1);
               S.drawMeshColor(S.tubeMesh(20), m.get(), S.gold);
            m.restore();
            //party hat rigid
            m.save();
               m.translate(0, .13* 2 + 0.05, 0);
               m.rotx(Math.PI/2);
               m.scale(0.05, 0.05, 0.05);
               S.drawMeshColor(S.coneMesh(100), m.get(), S.copper);
            m.restore();
         
            m.restore();

      m.restore();
   m.restore();

m.restore();
`,
events: `
   ;
`
};

}

