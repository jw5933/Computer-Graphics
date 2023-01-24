
rooms.raytrace5 = function() {

   lib3D();
   
   description = `Raytrace to quadrics<br>in a fragment shader
   <small>
      <p> 
      <input type=range id=red   value= 5> bg red
      <br> <input type=range id=green value=10> bg green
      <br> <input type=range id=blue  value=50> bg blue
      <br> <input type=range id=refract value=33> refract1
      <div id=iorInfo>&nbsp;</div>
      </p>
      <p>
      <b>Fog</b>
      <br> <input type=range id=fogRed   min= 50 value= 50> red
      <br> <input type=range id=fogGreen min = 50 value=80> green
      <br> <input type=range id=fogBlue  min = 50 value=100> blue
      <br> <input type=range id=fogDensity value= 0.5 min = 0 max = 1 step = 0.01> density
      </p>
   </small>
   `;
   
   code = {
   'init':`
   
      // DEFINE MATERIALS TO BE RENDERED VIA PHONG REFLECTANCE MODEL
   
      S.redPlastic    = [.2,.1,.1,0,  .5,.2,.2,0,  2,2,2,20,  0,0,0,0];
      S.greenPlastic  = [.1,.2,.1,0,  .2,.5,.2,0,  2,2,2,20,  0,0,0,0];
      S.bluePlastic   = [.1,.1,.2,0,  .2,.2,.5,0,  2,2,2,20,  0,0,0,0];
      S.whitePlastic  = [.2,.2,.2,0,  .5,.5,.5,0,  2,2,2,20,  0,0,0,0];
   `,
   
   fragment: `
   S.setFragmentShader(\`
   
      // DECLARE CONSTANTS, UNIFORMS, VARYING VARIABLES
      uniform float uTime;
      
      const int nQ = \` + S.nQ + \`;
      const int nL = \` + S.nL + \`;
      uniform vec3 uBgColor;
      uniform vec3 uLd[nL];
      uniform vec3 uLc[nL];
      uniform mat4 uQ[nQ];
      uniform mat4 uPhong[nQ];
      uniform int  uShape[nQ];
      uniform float uIor; //index of refraction
      
      uniform vec3 uFogColor;
      uniform float uFogDensity;
   
      varying vec3 vPos;
      float fl = 3.;
   
   /********* PSEUDO-CODE IMPLEMENTATION OF FRAGMENT SHADER **********/
   
   /*Compute surface normal: P, Q => N
   
      vec3 normalQ(vec3 P, mat4 Q)
   
         Just like in the course notes.
   */
   
      vec3 normalQ(vec3 P, mat4 Q) {
         /*
         a  b  c  d
         e  f  g  h
         i  j  k  l
         m  n  o  p
   
         a  e  i  m
         b  f  j  n
         c  g  k  o
         d  h  l  p
         */
         float x = P.x, y = P.y, z = P.z;
         float a = Q[0].x, e = Q[1].x, i = Q[2].x, m = Q[3].x;
         float b = Q[0].y, f = Q[1].y, j = Q[2].y, n = Q[3].y;
         float c = Q[0].z, g = Q[1].z, k = Q[2].z, o = Q[3].z;
         float d = Q[0].w, h = Q[1].w, l = Q[2].w, p = Q[3].w;
   
         
         float fx = 2.0*a*x + (b+e)*y + (c+i)*z + (d+m);
         float fy = (b+e)*x + 2.0*f*y + (g+j)*z + (h+n);
         float fz = (c+i)*x + (g+j)*y + 2.0*k*z + (l+o);
   
         return normalize(vec3(fx,fy,fz));
      }
   
   
   /*Trace a ray to a quadric: V, W, Q => [ tIn, tOut ]
   
      vec2 rayQ(vec3 V, vec3 W, mat4 Q)
   
         Just like in the course notes:
   
            First add homogeneous coordinates:
   
               V1 = vec4(V,1)
               W0 = vec4(W,0)
   
            Then compute quadratic equation:
   
               a: W0 . Q*W0
               b: V1 . Q*W0 + W0 . Q*V1
               c: V1 . Q*V1
   
            Then solve quadratic equation.
   
            Return both roots as a vec2
   */
   
      vec2 rayQ (vec4 V, vec4 W, mat4 Q){
         //W^T * Q * W
         float a = dot(W * Q, W);
         
         //V^T * Q * W  +  W^T * Q * V
         float b = dot(V * Q, W) + dot(W * Q, V);
         
         //V^T Q V
         float c = dot(V * Q, V);
         
         //solve quadratic
         //t = (-b ± √b^2 - 4ac) / 2a
         float disc = b * b - (4.0 * a * c);
   
         //return early if none (no real roots)
         if (disc < 0.0) return vec2(1000.0, -1000.0);
   
         float sq = sqrt(disc);
         //first point (which one should be the "first point"???? the smaller one)
         float t1 = (-b - sq) / (2.0 * a);
         //second point
         float t2 = (-b + sq) / (2.0 * a);
         
         return vec2(min(t1, t2), max(t1,t2));
      }
   
   
   /*
      Trace a ray to an intersection of quadric surfaces:
   
      Q1: T=[n,tIn,tOut], n, t.xy => T
      
         tIn  = t.x
         tOut = t.y
         if tIn > 0 and tIn < tOut and tIn < T.y
            T = [n,t]
         return T
      
      Q2: T=[n,tIn,tOut], n, t0.xy, t1.xy => T
      
         tIn  = max(t0.x,t1.x)
         tOut = min(t0.y,t1.y)
         if tIn > 0 and tIn < tOut and tIn < T.y
            i = t0.x==tIn ? 0 : 1
            T = [n+i,t]
         return T
      
      Q3: T=[n,tIn,tOut], n, t0.xy, t1.xy, t2.xy => T
      
         tIn  = max(t0.x,t1.x,t2.x)
         tOut = min(t0.y,t1.y,t2.y)
         if tIn > 0 and tIn < tOut and tIn < T.y
            i = t0.x==tIn ? 0 : t1.x==tIn ? 1 : 2
            T = [n+i,t]
         return T
      
      Q4: T=[n,tIn,tOut], n, t0.xy, t1.xy, t2.xy, t3.xy => T
      
         tIn  = max(t0.x,t1.x,t2.x,t3.x)
         tOut = min(t0.y,t1.y,t2.y,t3.y)
         if tIn > 0 and tIn < tOut and tIn < T.y
            i = t0.x==tIn ? 0 : t1.x==tIn ? 1 : t2.x==tIn ? 2 : 3
            T = [n+i,t]
         return T
   */
      //T = n, tIn, tOut
   
      vec3 rayToIntersect1(vec3 T, int n, vec2 t0){
         float tIn = t0.x;
         float tOut = t0.y;
   
         if (tIn > 0. && tIn < tOut && tIn < T.y){
            T = vec3(n, vec2(tIn, tOut));
         }
         return T;
      }
   
      vec3 rayToIntersect2(vec3 T, int n, vec2 t0, vec2 t1){
         float tIn = max(t0.x, t1.x);
         float tOut = min(t0.y, t1.y);
   
         if (tIn > 0. && tIn < tOut && tIn < T.y){
            int i = t0.x == tIn ? 0 : t1.x==tIn ? 1 : 2;
            T = vec3(n+i, vec2(tIn, tOut));
         }
         return T;
      }
   
      vec3 rayToIntersect3(vec3 T, int n, vec2 t0, vec2 t1, vec2 t2){
         float tIn  = max(max(t0.x, t1.x), t2.x);
         float tOut = min(min(t0.y, t1.y), t2.y);
   
         if (tIn > 0. && tIn < tOut && tIn < T.y){
            int i = t0.x == tIn ? 0 : t1.x == tIn ? 1 : 2;
            T = vec3(n+i, vec2(tIn, tOut));
         }  
         return T;
      }
   
      vec3 rayToIntersect4(vec3 T, int n, vec2 t0, vec2 t1, vec2 t2, vec2 t3){
         float tIn  = max(max(max(t0.x,t1.x),t2.x), t3.x);
         float tOut = min(min(min(t0.y,t1.y),t2.y), t3.y);
   
         if (tIn > 0. && tIn < tOut && tIn < T.y){
            int i = t0.x == tIn ? 0 : t1.x == tIn ? 1 : t2.x == tIn ? 2 : 3;
            T = vec3(n+i, vec2(tIn, tOut));
         }
         return T;
      }
   
   
      //vec2 rayQ (vec4 V, vec4 W, mat4 Q)
      vec3 rayScene(vec3 V, vec3 W){
         vec4 V1 = vec4(V, 1);
         vec4 W0 = vec4(W, 0);
         vec3 T = vec3(-1, 1000, 0);
   
         for(int n = 0; n < nQ; n++){
            if (uShape[n] == 0) continue;
            
            if (uShape[n] == 1)
               T = rayToIntersect1(T, n, rayQ(V1, W0, uQ[n]));
            else if (uShape[n] == 2)
               T = rayToIntersect2(T, n, rayQ(V1, W0, uQ[n]), rayQ(V1, W0, uQ[n+1]) );
            else if (uShape[n] == 3)
               T = rayToIntersect3(T, n, rayQ(V1, W0, uQ[n]), rayQ(V1, W0, uQ[n+1]), rayQ(V1, W0, uQ[n+2]) );
            else if (uShape[n] == 4)
               T = rayToIntersect4(T, n, rayQ(V1, W0, uQ[n]), rayQ(V1, W0, uQ[n+1]), rayQ(V1, W0, uQ[n+2]), rayQ(V1, W0, uQ[n+3]) );
         }
         return T;
      }
   
      vec3 shadeSurface(vec3 P, vec3 N, mat4 M){
         vec3 ambient  = M[0].rgb;
         vec3 diffuse  = M[1].rgb;
         vec3 specular = M[2].rgb;
         float p       = M[2].a;
   
         vec3 c = mix(ambient, uBgColor, .3);
         vec3 E = vec3(0.,0.,1.); //eye
   
         for (int l = 0 ; l < nL ; l++) {
            
            float t = -1.;
            
            // Shadows
            // for (int n = 0 ; n < nS ; n++)
            // t = max(t, raySphere(P, uLd[l], uS[n]));
   
            // Lighting if not in shadow
            if (t < 0.) {
               vec3 R = 2. * dot(N, uLd[l]) * N - uLd[l];
               c += uLc[l] * (diffuse * max(0.,dot(N, uLd[l]))
                           + specular * pow(max(0., dot(R, E)), p));
            }
         }
   
         return c;
      }
   
      vec3 refractRay(vec3 W, vec3 N, float n1){
         vec3 C1 = N * dot(W, N);
         vec3 S1 = W - C1;
   
         float sin1 = length(S1); 
         float cos1 = sqrt(1.0 - sin1 * sin1); // incoming angle
         
         float sin2 = sin1 * n1;
         float cos2 = sqrt(1.0 - sin2 * sin2); //outgoing angle
   
         //components of ray
         vec3 C2 = C1 * cos2 / cos1;
         vec3 S2 = S1 * sin2 / sin1;
   
         //direction of outgoing ray
         vec3 W2 = C2 + S2;
         return W2;
      }
   
      void main() {
         vec3 color = uBgColor;
   
         // Main loop
         vec3 V = vec3(0.,0., fl);
         vec3 W = normalize(vec3(vPos.xy, -fl));
         vec3 T = rayScene(V, W);
         vec3 firstT = T;
   
         mat4 Q;
         mat4 phong1;
         mat4 phong2;
         mat4 phong3;
         
         int n = int(T.x); // REMEMBER, YOU NEED TO MAKE A LOOP HERE, AS SHOWN ABOVE.
         
         for (int i = 0; i < nQ; i++){
            if (i == n){
               Q = uQ[i];
               phong1 = uPhong[i];
            }
         }
   
         if (n >= 0){
            //    Compute surface point P = V + T.y * W
            vec3 P = V + T.y * W;
            //    Shade with Phong
            vec3 N = normalQ(P, Q);
            color = shadeSurface(P, N, phong1);
   
            //    Do reflection:
            vec3 R = 2. * dot(N, -W) * N + W;
            T = rayScene(P, R);
            n = int(T.x);
            for (int i = 0; i < nQ; i++){
               if (i == n){
                  Q = uQ[i];
                  phong1 = uPhong[i];
               }
            }
            
            if (n >= 0){
               vec3 M = P + T.y * W;
               color += shadeSurface(M, normalQ(M, Q), phong1)/3.;
            }
   
   
            // Do refraction:
            // (1) SHOOT RAY TO FIND REAR OF THIS OBJECT (USE 2nd ROOT):
            W = refractRay(W, N, uIor);
            T = rayScene (P - .001*W, W); // [n,tIn,tOut]
            
            n = int(T.x);
            for (int i = 0; i < nQ; i++){
               if (i == n){
                  Q = uQ[i];
                  phong2 = uPhong[i];
               }
            }
   
            //find outside objects
            P = P + T.z * W; // point at rear (T.z is tOut)
            N = normalQ(P, Q);
   
            // (2) SHOOT RAY FROM REAR OF THIS OBJECT TO FIND ANOTHER OBJECT:
            W = refractRay (W, N, 1.0/uIor);
            T = rayScene(P , W);
            
            n = int(T.x);
            for (int i = 0; i < nQ; i++){
               if (i == n){
                  Q = uQ[i];
                  phong3 = uPhong[i];
               }
            }
   
            if (n >= 0){
               vec3 M = P + T.y * W; // point on other surface (T.y is tIn) 
               color += phong2[1].rgb *
                        shadeSurface(M, normalQ(M, Q), phong3);
            }
         }
         
         //add fog
         // float f = pow(uFogDensity, firstT.y);
         // color = mix(uFogColor, color, f);
         
         //add fog
         vec3 p = vPos + vec3(.1* uTime *2.0, 0., .1*uTime);
         p*= 1.5;
         float fn = fractal(p);
         color = mix(mix((uFogColor - 0.2), uFogColor, fn), color, pow(uFogDensity,  firstT.y + 1.5 * fn));
      
         gl_FragColor = vec4(sqrt(color), 1.);
      }
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
   render: `
   
      // USEFUL VECTOR FUNCTIONS
   
      let add = (a,b) => [ a[0]+b[0], a[1]+b[1], a[2]+b[2] ];
      let dot = (a,b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
      let norm = v => Math.sqrt(dot(v,v));
      let normalize = v => { let s = norm(v); return [ v[0]/s, v[1]/s, v[2]/s ]; }
      let scale = (v,s) => [ s * v[0], s * v[1], s * v[2] ];
      let subtract = (a,b) => [ a[0]-b[0], a[1]-b[1], a[2]-b[2] ];
   
      // SEND LIGHT SOURCE DATA TO GPU
   
      let ldData = [ normalize([1,1,1]),
                     normalize([-1,-1,-1]) ];
      S.setUniform('3fv', 'uLd', ldData.flat());
      S.setUniform('3fv', 'uLc', [ 1,1,1, .5,.3,.1 ]);
   
      // DEFINE NUMBER OF LIGHTS FOR GPU
   
      S.nL = ldData.length;
   
      // SEND BACKGROUND COLOR TO GPU
   
      S.setUniform('3fv', 'uBgColor', [ red.value   / 100,
                                        green.value / 100,
                                        blue.value  / 100 ]);
   
      // SEND INDEX OF REFRACTION TO GPU
   
      let ior = refract.value / 100 + 1;
      S.setUniform('1f', 'uIor', ior);
      S.setUniform('1f', 'uTime', time);
   
      // SEND FOG TO GPU
      S.setUniform('3fv', 'uFogColor', [fogRed.value  /100,
                                       fogGreen.value/100,
                                       fogBlue.value /100 ]);
      S.setUniform('1f', 'uFogDensity', fogDensity.value);
   
      // DIFFERENT QUADRIC SURFACES
   
   //                xx        yy         zz           c
   
      let qSlabX  = [1,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,-1]; // x*x - 1 <= 0
      let qSlabY  = [0,0,0,0, 0,1,0,0, 0,0,0,0, 0,0,0,-1]; // y*y - 1 <= 0
      let qSlabZ  = [0,0,0,0, 0,0,0,0, 0,0,1,0, 0,0,0,-1]; // z*z - 1 <= 0
      let qSphere = [1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,-1]; // x*x + y*y + z*z - 1 <= 0
      let qTubeX  = [0,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,-1]; // y*y + z*z - 1 <= 0
      let qTubeY  = [1,0,0,0, 0,0,0,0, 0,0,1,0, 0,0,0,-1]; // x*x + z*z - 1 <= 0
      let qTubeZ  = [1,0,0,0, 0,1,0,0, 0,0,0,0, 0,0,0,-1]; // x*x + y*y - 1 <= 0
      let qParaboloid = [1, 0, 0, 0, 	//c1
      0, 1, 0, 0, 	//c2
      0, 0, 0, 0, 	//c3
      0, 0, -1, 0];	//c4
      let qHyperboloid = //hyperboloid of one sheet
      [1, 0, 0, 0,
      0, 1, 0, 0,
      0, 0, -1, 0,
      0, 0, 0, -1];
      let qCone = //cone
      [1, 0, 0, 0,
      0, 1, 0, 0,
      0, 0, -1, 0,
      0, 0, 0, 0];
      let qHyperboloid2 = //hyperboloid of two sheets
      [1, 0, 0, 0,
      0, 1, 0, 0,
      0, 0, -1, 0,
      0, 0, 0, 1];
   
      // SHAPES ARE INTERSECTIONS OF QUADRIC SURFACES
   
      let shape = [], coefs = [], xform = [], phong = [], M;
   
      let sphere = (m, M) => {
         shape.push(1);
         phong.push(m);
         xform.push(M);
         coefs.push(qSphere);
      }
   
      let tubeX = (m, M) => {
         shape.push(2, 0);
         phong.push(m, m);
         xform.push(M, M);
         coefs.push(qTubeX, qSlabX);
      }
   
      let tubeY = (m, M) => {
         shape.push(2, 0);
         phong.push(m, m);
         xform.push(M, M);
         coefs.push(qTubeY, qSlabY);
      }
   
      let tubeZ = (m, M) => {
         shape.push(2, 0);
         phong.push(m, m);
         xform.push(M, M);
         coefs.push(qTubeZ, qSlabZ);
      }
   
      let cube = (m, M) => {
         shape.push(3, 0, 0);
         phong.push(m, m, m);
         xform.push(M, M, M);
         coefs.push(qSlabX, qSlabY, qSlabZ);
      }
   
      let rand = (m, M1, M2, M3, M4) => {
         shape.push(4, 0, 0, 0);
         phong.push(m, m, m, m);
         xform.push(M1, M2, M3, M4);
         coefs.push(qSlabX, qSlabY, qSlabZ, qTubeX);
      }
   
      let octahedron = (m, M) => {
         shape.push(4, 0, 0, 0);
         phong.push(m, m, m, m);
         xform.push(M, M, M, M);
         coefs.push([1, 2, 2, 0,  0, 1, 2, 0,  0,0,1,0,  0,0,0,-1]);
         coefs.push([1,-2,-2, 0,  0, 1, 2, 0,  0,0,1,0,  0,0,0,-1]);
         coefs.push([1,-2, 2, 0,  0, 1,-2, 0,  0,0,1,0,  0,0,0,-1]);
         coefs.push([1, 2,-2, 0,  0, 1,-2, 0,  0,0,1,0,  0,0,0,-1]);
      }
   
      let hourglass = (m, M) => {
         shape.push(2, 0);
         phong.push(m, m);
         xform.push(M, M);
         coefs.push(qHyperboloid, qSlabZ);
      }
   
      // CREATE THE SCENE
   
      cube(S.bluePlastic,
            mScale(.1,.1,.1,
            mRoty(time * 1.1,
            mRotz(time * 1.2,
            matrixTranslate(-Math.sin(time)*.5,0,Math.cos(time)*.5+.5)))));
   
      cube(S.whitePlastic,
          mScale(.18,.18,.18,
          mRoty(time * 1.2,
          mRotx(time * 1.1,
          matrixTranslate(Math.sin(time)*.5, 0, -Math.cos(time)*2.)))));
   
      // octahedron(S.greenPlastic,
      //    mScale(.18,.18,.18,
      //    mRoty(1.2,
      //    mRotz(1.1,
      //    mRotx(1.3,
      //    matrixTranslate(0, 0, 0))))));
   
      // cube(S.whitePlastic,
      //    mScale(.18,.18,.18,
      //    matrixTranslate(0, 0, 0)));
      
      // cube(S.whitePlastic,
      //    mScale(.05,.3,.05,
      //    mRoty(1.2,
      //    mRotx(Math.PI/4 * Math.cos(time),
      //    mRotz(Math.PI,
      //    matrixTranslate(0, 0, 0.5))))));
   
         // cube(S.whitePlastic,
      //     mScale(.18,.18,.18,
      //     mRoty(time * 1.2,
      //     mRotz(time * 1.1,
      //     mRotx(time * 1.3,
      //     matrixTranslate(0,Math.cos(time)*.2,.5))))));
   
      // sphere(S.bluePlastic,
      //        mScale(.2,.15,.18,
      //        mRoty(time * 1.3,
      //        mRotz(time * 1.1,
      //        mRotx(time * 1.2,
      //        matrixTranslate(Math.sin(time)*.5,0,-Math.cos(time)*.5+.5))))));
   
      // SEND SCENE DATA TO GPU
   
      for (let n = 0 ; n < coefs.length ; n++) {
         let IM = matrixInverse(xform[n]);
         coefs[n] = matrixMultiply(matrixTranspose(IM), matrixMultiply(coefs[n], IM));
      }
   
      S.setUniform('1iv', 'uShape', shape);
      S.setUniform('Matrix4fv', 'uQ', false, coefs.flat());
      S.setUniform('Matrix4fv', 'uPhong', false, phong.flat());
   
      // DEFINE NUMBER OF QUADRIC SURFACES FOR GPU
   
      S.nQ = coefs.length;
   
      // RENDER THIS ANIMATION FRAME
   
      S.gl.drawArrays(S.gl.TRIANGLE_STRIP, 0, 4);
   
      // SET ANY HTML INFO
   
      iorInfo.innerHTML = 'index of refraction = ' + (ior * 100 >> 0) / 100;
   `,
   events: `
      ;
   `
   };
   
   }
   
   
   