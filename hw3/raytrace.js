
rooms.raytrace = function() {

    lib3D();
    
    description = `Raytracing to spheres<br>in a fragment shader
    <small>
        <p>
        <b>Background color</b>
        <br> <input type=range id=red   value= 5> red
        <br> <input type=range id=green value=10> green
        <br> <input type=range id=blue  value=50> blue
        </p>
        <p>
        <b>Fog</b>
        <br> <input type=range id=fogRed   min= 50 value= 50> red
        <br> <input type=range id=fogGreen min = 50 value=80> green
        <br> <input type=range id=fogBlue  min = 50 value=100> blue
        <br> <input type=range id=fogDensity value= 0.5 min = 0 max = 0.9 step = 0.01> density
        </p>
    </small>
    `;
    
    code = {
    'init':`
      // DEFINE NUMBER OF SPHERES AND NUMBER OF LIGHTS
   
      S.nS = 100;
      S.nP = 5;
      S.nL = 2;
      S.nQ = 2;

      //define geneneral 2nd order shapes
      S.defaultShapes = [
         //slab
         [0, 0, 0, 0,
         0, 0, 0, 0,
         0, 0, 1, 0,
         0, 0, 0, -1],
         //tube
         [1, 0, 0, 0,
         0, 1, 0, 0,
         0, 0, 0, 0,
         0, 0, 0, -1],
         //paraboloid
         [1, 0, 0, 0, 	//c1
         0, 1, 0, 0, 	//c2
         0, 0, 0, 0, 	//c3
         0, 0, -1, 0],	//c4
         //hyperboloid of one sheet
         [1, 0, 0, 0,
         0, 1, 0, 0,
         0, 0, -1, 0,
         0, 0, 0, -1],
         //cone
         [1, 0, 0, 0,
         0, 1, 0, 0,
         0, 0, -1, 0,
         0, 0, 0, 0],
         //hyperboloid of two sheets
         [1, 0, 0, 0,
         0, 1, 0, 0,
         0, 0, -1, 0,
         0, 0, 0, 1],
      ];
   
      // DEFINE MATERIALS TO BE RENDERED VIA PHONG REFLECTANCE MODEL
   
      let materials = [
         [.15,.05,.025,0, .3,.1,.05,0, .6,.2,.1,3, 0,0,0,0], // COPPER
         [.25,.15,.025,0, .5,.3,.05,0, 1,.6,.1,6,  0,0,0,0], // GOLD
         [.25,0,0,0,      .5,0,0,0,    2,2,2,20,   0,0,0,0], // PLASTIC
         [.05,.05,.05,0,  .1,.1,.1,0,  1,1,1,5,    0,0,0,0], // LEAD
         [.1,.1,.1,0,     .1,.1,.1,0,  1,1,1,5,    0,0,0,0], // SILVER
      ];
   
      S.material = [];
      for (let n = 0 ; n < S.nS ; n++)
         S.material.push(materials[n % materials.length]);
   
      // INITIALIZE SPHERE POSITIONS AND VELOCITIES
   
      S.sPos = [];
      S.sVel = [];
      for (let n = 0 ; n < S.nS ; n++) {
         S.sPos.push([ Math.random() - .5,
                     Math.random() - .5,
                     Math.random() - .5 ]);
         S.sVel.push([0,0,0]);
      }
    `,
    fragment: `
    S.setFragmentShader(\`
    
       // DECLARE CONSTANTS, UNIFORMS, VARYING VARIABLES
    
       const int nS = \` + S.nS + \`;
       const int nL = \` + S.nL + \`;
       const int nP = \` + S.nP + \`;
       const int nQ = \` + S.nQ + \`;

       const int polyCount = 5;
       const mat4 matrixIdentity = mat4(
        1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1, 0,
        0, 0, 0, 1
       );

       uniform float uTime;
       uniform vec3 uBgColor;
       uniform vec3 uFogColor;
       uniform vec2 uFogDensity;

       uniform mat4 uShapes[3];
       uniform mat4 uQ[nQ];
       uniform mat4 uQM;
       uniform mat4 uIM[3];
       
       uniform vec4 uS[nS];
       uniform mat4 uMaterial[nS];
       uniform vec4 uCube[6];
       uniform vec4 uOctahedron[8];
       uniform vec4 uPyramid[5];
       uniform vec3 uLd[nL];
       uniform vec3 uLc[nL];
    
       varying vec3 vPos;
    
       // DEFINE CAMERA FOCAL LENGTH
    
       float fl = 2.;
      
      mat4 transpose(mat4 m) {
         vec4 c0 = m[0];
         vec4 c1 = m[1];
         vec4 c2 = m[2];
         vec4 c3 = m[3];
     
         mat4 transposed = mat4(
                     vec4(c0.x, c1.x, c2.x, c3.x),
                     vec4(c0.y, c1.y, c2.y, c3.y),
                     vec4(c0.z, c1.z, c2.z, c3.z),
                     vec4(c0.w, c1.w, c2.w, c3.w)
                     );
         return transposed;
      }

      float turbulence(vec3 p) {
         float t = 0., f = 1.;
         for (int i = 0 ; i < 4 ; i++) {
            t += abs(noise(f * p)) / f;
            f *= 2.;
         }
         return t;
      }

      
       //trace to general 2nd order
      vec2 rayShape (vec4 V, vec4 W, mat4 Q){
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

      vec3 quadricSurfaceNormal(mat4 Q, vec3 P) {
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
       
       // TRACE A RAY TO A HALFSPACE
       float rayHalfspace(vec3 V, vec3 W, vec4 H) {
          vec4 V1 = vec4(V, 1.);
          vec4 W0 = vec4(W, 0.);
          return -dot(H, V1) / dot(H, W0);
       }
    

       // TRACE A RAY TO A TRANSFORMED CUBE
        vec4 rayCube(vec3 V, vec3 W, mat4 IM) {
            vec3 N = vec3(0.);
            float tIn = -1000., tOut = 1000.;
            for (int i = 0 ; i < 6 ; i++) {
                vec4 H = uCube[i] * IM;
            H /= sqrt(dot(H.xyz, H.xyz));
            float t = rayHalfspace(V, W, H);
            if (dot(W, H.xyz) < 0.) {
                if (t > tIn)
                N = H.xyz;
                tIn = max(tIn, t);
            }
            else
                tOut = min(tOut, t);
            }
            return vec4(N, tIn < tOut ? tIn : -1.);
        }
    
       //trace ray to a transformed octahedron
        vec4 rayOctahedron(vec3 V, vec3 W, mat4 IM) {
            vec3 N = vec3(0.);
            float tIn = -1000., tOut = 1000.;
            
            for (int i = 0 ; i < 8 ; i++) {
                vec4 H = uOctahedron[i] * IM;
                H /= sqrt(dot(H.xyz, H.xyz));
                float t = rayHalfspace(V, W, H);
                
                if (dot(W, H.xyz) < 0.) {
                    if (t > tIn)
                    N = H.xyz;
                    tIn = max(tIn, t);
                }
                else
                    tOut = min(tOut, t);
                }
                return vec4(N, tIn < tOut ? tIn : -1.);
            }

         //trace ray to pyramid
         vec4 rayPyramid(vec3 V, vec3 W, mat4 IM) {
            vec3 N = vec3(0.);
            float tIn = -1000., tOut = 1000.;
            for (int i = 0; i < 5; i++) {
               vec4 H = uPyramid[i] * IM;
               H /= sqrt(dot(H.xyz, H.xyz));
               float t = rayHalfspace(V, W, H);
               if (dot(W, H.xyz) < 0.) {
                  if (t > tIn)
                  N = H.xyz; //normal
                  tIn = max(tIn, t); //last t in
               }
               else
                  tOut = min(tOut, t); //first t out
               }
               return vec4(N, tIn < tOut ? tIn : -1.); //if t in > t out then we've hit poly
        }

       // TRACE A RAY TO A SPHERE
       float raySphere(vec3 V, vec3 W, vec4 S) {
         V -= S.xyz;
         V += .01 * W;
         float b = dot(V, W);
         float d = b * b - dot(V, V) + S.w * S.w;
         return d < 0. ? -1. : -b - sqrt(d);
      }

      // ==================== SHADING ========================
      // SHADE A SHAPE
      vec3 shadeSurface(vec3 P, vec3 W, vec3 N, mat4 M) {
         vec3 ambient  = M[0].rgb;
         vec3 diffuse  = M[1].rgb;
         vec3 specular = M[2].rgb;
         float p       = M[2].a;
   
         vec3 c = mix(ambient, uBgColor, .1);
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
    
      // SHADE A SPHERE AT ONE SURFACE POINT
      vec3 shadeSphere(vec3 P, vec4 S, mat4 M) {
   
         // EXTRACT PHONG PARAMETERS FROM MATERIAL MATRIX
   
         vec3  ambient  = M[0].rgb;
         vec3  diffuse  = M[1].rgb;
         vec3  specular = M[2].rgb;
         float p        = M[2].a;
   
         // COMPUTE NORMAL, INIT COLOR, APPROXIMATE VECTOR TO EYE
   
         vec3 N = normalize(P - S.xyz);
         vec3 c = mix(ambient, uBgColor, .3);
         vec3 E = vec3(0.,0.,1.);
   
         // LOOP THROUGH LIGHT SOURCES
   
         for (int l = 0 ; l < nL ; l++) {

            // ADD DIFFUSE AND SPECULAR ONLY IF NOT IN SHADOW

            float t = -1.;
            for (int n = 0 ; n < nS ; n++){
               t = max(t, raySphere(P, uLd[l], uS[n]));
            }
            for (int n = 0 ; n < nP ; n++){
               t = max(t, rayOctahedron(P, uLd[l], uIM[n]).w);
            }

            // COMPUTE DIFFUSE AND SPECULAR FOR THIS LIGHT SOURCE

            if (t < 0.) {
               vec3 R = 2. * dot(N, uLd[l]) * N - uLd[l];
               c += uLc[l] * (diffuse * max(0.,dot(N, uLd[l]))
                           + specular * pow(max(0., dot(R, E)), p));
            }
         }
         c *= 1. + .5 * noise(3.*N); // OPTIONAL SPOTTY TEXTURE
         return c;
      }

      // SHADE A POLY AT ONE SURFACE POINT
      vec3 shade(vec3 P, vec3 W, vec3 N, mat4 M) {
         // EXTRACT PHONG PARAMETERS FROM MATERIAL MATRIX
         vec3  ambient  = M[0].rgb;
         vec3  diffuse  = M[1].rgb;
         vec3  specular = M[2].rgb;
         float p        = M[2].a;
   
         // COMPUTE NORMAL, INIT COLOR, APPROXIMATE VECTOR TO EYE
         vec3 c = mix(ambient, uBgColor, .3); //colour
         vec3 E = vec3(0.,0.,1.); //?
   
         // LOOP THROUGH LIGHT SOURCES
         for (int l = 0 ; l < nL ; l++) {
         // ADD DIFFUSE AND SPECULAR ONLY IF NOT IN SHADOW
         
         float t = -1.;
         for (int n = 0 ; n < nS ; n++)
               t = max(t, raySphere(P, uLd[l], uS[n]));
         

         // COMPUTE DIFFUSE AND SPECULAR FOR THIS LIGHT SOURCE where its not in shadow
         if (t < 0.) {
            vec3 R = 2. * dot(N, uLd[l]) * N - uLd[l];
            vec3 reflectDir = reflect(-uLd[l], N);
            c += uLc[l] * (
               diffuse * max(0.,dot(N, uLd[l]))
               + specular * pow(max(dot(E, reflectDir), 0.0), p)
               );
         }
      }
      return c;
      }

    
      void main() {
   
         // BACKGROUND COLOR IS THE DEFAULT COLOR
   
         vec3 color = uBgColor;
   
         // DEFINE RAY INTO SCENE FOR THIS PIXEL
   
         vec3 V = vec3(0.,0.,fl);
         vec3 W = normalize(vec3(vPos.xy, -fl));
   
         
         float tMin = 10000.;
         mat4 Q;
         vec2 t = vec2(-1000., 1000.0);
         
         // loop through surfaces
         for (int n = 0; n < nQ; n++) {
            vec2 ts = rayShape(vec4(V, 1), vec4(W,0), uQ[n]);
            
            /*
            if the first point of entry is greater than prev surface, 
            we want to draw this one
            */
           
            if (ts.x > t.x){
               Q = uQ[n];

               //intersection
               t.x = max(ts.x, t.x);
               t.y = min(ts.y, t.y);
            }

            //draw if first object hit, and draw only surface from exterior
            if (t.x < tMin && t.x < t.y) {
               vec3 P = V + t.x * W;
               vec3 N = quadricSurfaceNormal(Q, P);
               color = shadeSurface(P, W, N, uMaterial[0]);
               tMin = t.x;
            }
            else{
               vec4 Nt = rayOctahedron(V, W, uIM[0]);
               if (Nt.w > 0. && Nt.w < tMin) {
                  vec3 a = mix(vec3(.1), uBgColor, .3);
                  color = a + vec3(max(0., dot(Nt.xyz, vec3(.5))));
               }
            }
         }

         // LOOP THROUGH SPHERES
         for (int n = 0 ; n < nS ; n++) {
            float t = raySphere(V, W, uS[n]);
            if (t > 0. && t < tMin) {
               // IF THIS IS THE NEAREST SPHERE, DO PHONG SHADING
   
               vec3 P = V + t * W;
               color = shadeSphere(P, uS[n], uMaterial[n]);
               tMin = t;
   
               // SEE WHETHER ANY OTHER SPHERE IS VISIBLE VIA REFLECTION
   
               vec3 N = normalize(P - uS[n].xyz);
               vec3 R = 2. * dot(N, -W) * N + W;
               float rtMin = 10000.;
               vec3 rColor;
               for (int rn = 0 ; rn < nS ; rn++) {
                     float rt = raySphere(P, R, uS[rn]);
                     if (rt > 0. && rt < rtMin) {
                        rtMin = rt;
                        rColor = shadeSphere(P + rt * R, uS[rn], uMaterial[rn]);
                     }
               }
               if (rtMin < 10000.)
                  color += .5 * rColor;
            }
         }
         
         // RAY TRACE TO THE POLY
         // for (int n = 0; n < nP; n++){
         //    vec4 Nt = rayOctahedron(V, W, uIM[n]);
         //    vec3 P = V + Nt.w * W;
            
         //    if (Nt.w > 0. && Nt.w < tMin) {
         //    //shading poly
         //    color = shade(P, W, Nt.xyz, uMaterial[1]);

         //    //SEE WHETHER ANY OTHER SPHERE IS VISIBLE VIA REFLECTION
         //    for (int sn = 0; sn < 8; sn ++){
         //       vec3 N = normalize(P - uOctahedron[sn].xyz);
         //       vec3 R = 2. * dot(N, -W) * N + W;
         //       float rtMin = 10000.;
         //       vec3 rColor;

         //       //reflect spheres
         //       for (int rn = 0 ; rn < nS ; rn++) {
         //          float rt = raySphere(P, R, uS[rn]);
         //          if (rt > 0. && rt < rtMin) {
         //             rtMin = rt;
         //             rColor = shadeSphere(P + rt * R, uS[rn], uMaterial[rn]);
         //          }
         //       }
         //       //reflect
         //       if (rtMin < 10000.)
         //          color += .25 * rColor;
         //       }
         //       //color *= 1. + .5 * noise(3.* N); // OPTIONAL SPOTTY TEXTURE
         //    }
         // }
         
         //add fog
         vec3 p = vPos + vec3(.1* uTime *2.0, 0., .1*uTime);
         p*= 1.5;
         float fn = fractal(p);
         color = mix(mix((uFogColor - 0.2), uFogColor, fn), color, pow(uFogDensity.x,  tMin + 1.5 * fn));
   
         // SET PIXEL COLOR
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
    
       // HANDY DANDY VECTOR LIBRARY
    
       let add = (a,b) => [ a[0]+b[0], a[1]+b[1], a[2]+b[2] ];
       let dot = (a,b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
       let norm = v => Math.sqrt(dot(v,v));
       let normalize = v => {
          let s = Math.sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
          return [ v[0]/s, v[1]/s, v[2]/s ];
       }
       let scale = (v,s) => [ s * v[0], s * v[1], s * v[2] ];
       let subtract = (a,b) => [ a[0]-b[0], a[1]-b[1], a[2]-b[2] ];
    
       // DEFINE RADIUS OF EACH SPHERE
    
       let radius = .07;
    
       // SEND LIGHT SOURCE DATA TO GPU
    
       let ldData = [ normalize([1,1,1]),
                      normalize([-1,-1,-1]) ];
       S.setUniform('3fv', 'uLd', ldData.flat());
       S.setUniform('3fv', 'uLc', [ 1,1,1, .5,.3,.1 ]);
    
       // SEND ANIMATION TIME TO GPU
    
       S.setUniform('1f', 'uTime', time);
    
       // MOVE SPHERES INTO PLACE FOR THIS ANIMATION FRAME
    
       for (let n = 0 ; n < S.nS ; n++) {
          for (let i = 0 ; i < 3 ; i++) {
             S.sVel[n][i] += .003 * Math.cos(time + (2+i) * n);
             S.sVel[n][i] += .01 * (Math.random() - .5);
             S.sPos[n][i] += .1 * S.sVel[n][i];
          }
          S.sPos[n] = scale(normalize(S.sPos[n]), .7);
       }
    
       // AVOID SPHERE INTERPENETRATION
    
       for (let m = 0 ; m < S.nS ; m++)
       for (let n = 0 ; n < S.nS ; n++)
          if (m != n) {
             let D = subtract(S.sPos[m], S.sPos[n]);
             let d = norm(D);
             if (d < 2 * radius) {
                let t = 2 * radius - d;
                for (let i = 0 ; i < 3 ; i++) {
                   S.sPos[m][i] += t * D[i] / d;
                   S.sPos[n][i] -= t * D[i] / d;
                }
             }
          }
    
       // SEND SPHERES DATA TO GPU
    
       let sData = [];
       for (let n = 0 ; n < S.nS ; n++)
          sData.push(S.sPos[n], radius);
       S.setUniform('4fv', 'uS', sData.flat());
       S.setUniform('Matrix4fv', 'uMaterial', false, S.material.flat());
    
       // SEND BACKGROUND COLOR TO GPU
    
       S.setUniform('3fv', 'uBgColor', [ red.value  /100,
                                         green.value/100,
                                         blue.value /100 ]);
      
       // SEND FOG TO GPU
       S.setUniform('3fv', 'uFogColor', [ fogRed.value  /100,
                                         fogGreen.value/100,
                                         fogBlue.value /100 ]);
      S.setUniform('2fv', 'uFogDensity', [fogDensity.value, fogDensity.value]);
      // console.log(fogDensity.value);
      
      // BEGINNINGS OF A MATRIX LIBRARY
    
      let matrixIdentity = () => 
         [1,0,0,0, 
         0,1,0,0, 
         0,0,1,0, 
         0,0,0,1];

      let matrixTranslate = (x,y,z) => 
         [1,0,0,0, 
         0,1,0,0, 
         0,0,1,0, 
         x,y,z,1]; //tube -> down is right -> left

      let matrixScale = (x,y,z) => 
         [x,0,0,0, 
         0,y,0,0, 
         0,0,z,0, 
         0,0,0,1];

      let matrixRotx = t => {
         let c = Math.cos(t), s = Math.sin(t);
         return [1,0,0,0, 0,c,s,0, 0,-s,c,0, 0,0,0,1];
      }
      let matrixRoty = t => {
         let c = Math.cos(t), s = Math.sin(t);
         return [c,0,-s,0, 0,1,0,0, s,0,c,0, 0,0,0,1];
      }
      let matrixRotz = t => {
         let c = Math.cos(t), s = Math.sin(t);
         return [c,s,0,0, -s,c,0,0, 0,0,1,0, 0,0,0,1];
      }
   
      // RENDER THE POLYHEDRON
      let polyM = [];
      let iM = matrixIdentity();

      let anim = iM;
      anim = matrixMultiply(anim, matrixScale(0.1, 0.1, 0.1));
      
      polyM.push(matrixInverse(anim));
      S.setUniform('Matrix4fv', 'uIM', false, polyM.flat());

      //shapes
      let qIMs = [];
      let qM = matrixIdentity();
      // qM = matrixMultiply(qM, matrixTranslate(.5,0,0));
      qM = matrixMultiply(qM, matrixRotx(time));
      qM = matrixMultiply(qM, matrixRoty(Math.PI/2));
      qM = matrixMultiply(qM, matrixRotx(2*time));
      qM = matrixMultiply(qM, matrixScale(.1,.1,.1));
      let qIM = matrixInverse(qM);

      qIMs.push(qIM);
      qIMs.push(qIM);

      qM = matrixIdentity();
      // qM = matrixMultiply(qM, matrixRotx(time));
      // qM = matrixMultiply(qM, matrixRotz(time));
      // qM = matrixMultiply(qM, matrixRotx(2*time));
      qM = matrixMultiply(qM, matrixScale(.2,.2,.1));
      
      // qIM = matrixInverse(qM);
      // qIMs.push(qIM);

      let wQs = [1, 0, 1];
      let qs = [];
      for (let n = 0 ; n < S.nQ; n++) {
         let Q = S.defaultShapes[wQs[n]];
         qs = qs.concat(
            matrixMultiply(matrixTranspose(qIMs[n]),
                           matrixMultiply(Q, qIMs[n])));
      }
      S.setUniform('Matrix4fv', 'uQ', false, qs);
      
      S.setUniform('Matrix4fv', 'uQM', false, qM);


      S.setUniform('4fv', 'uCube', [
         -1,0,0,-1, 1,0,0,-1,
         0,-1,0,-1, 0,1,0,-1,
         0,0,-1,-1, 0,0,1,-1,
      ]);
        
      S.setUniform('4fv', 'uOctahedron', [
         -1,-1,-1,-1, 1,1,1,-1,
         1,-1,-1,-1, -1,1,1,-1,
         -1,1,-1,-1, 1,-1,1,-1,
         -1,-1,1,-1, 1,1,-1,-1,
      ]);

      S.setUniform('4fv', 'uPyramid', [
          1, 1,  1, -1,
          1, 1, -1, -1,
         -1, 1,  1, -1,
         -1, 1, -1, -1,
          0,-1,  0, -1
      ]);
    
      // RENDER THIS ANIMATION FRAME
   
      S.gl.drawArrays(S.gl.TRIANGLE_STRIP, 0, 4);
    `,
    events: `
       ;
    `
    };
    
    }
    
    
    