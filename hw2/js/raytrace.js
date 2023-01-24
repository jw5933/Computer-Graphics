
rooms.raytrace = function() {

    lib3D();
    
    description = `Raytracing to spheres<br>in a fragment shader
    <small>
        <p>
        <b>Background color</b>
        <br> <input type=range id=red   value= 5> red
        <br> <input type=range id=green value=10> green
        <br> <input type=range id=blue  value=50> blue
    </small>
    `;
    
    code = {
    'init':`
    
       // DEFINE NUMBER OF SPHERES AND NUMBER OF LIGHTS
    
       S.nS = 20;
       S.nP = 5;
       S.nL = 2;
    
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
       const int polyCount = 5;
       const mat4 matrixIdentity = mat4(
        1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1, 0,
        0, 0, 0, 1
       );
       uniform float uTime;
       uniform vec3 uBgColor;
       uniform vec4 uS[nS];
       uniform mat4 uSm[nS];
       uniform mat4 uIM[nP];
       uniform vec4 uCube[6];
       uniform vec4 uOctahedron[8];
       uniform vec4 uPyramid[5];
       uniform vec3 uLd[nL];
       uniform vec3 uLc[nL];
    
       varying vec3 vPos;
    
       // DEFINE CAMERA FOCAL LENGTH
    
       float fl = 3.;
    
       
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


      //SHADE A POLY AT ONE SURFACE POINT
      vec3 shadePoly(vec3 P, vec4 S, mat4 M, vec3 norm) {
         // EXTRACT PHONG PARAMETERS FROM MATERIAL MATRIX
   
         vec3  ambient  = M[0].rgb;
         vec3  diffuse  = M[1].rgb;
         vec3  specular = M[2].rgb;
         float p        = M[2].a*2.0;
   
         vec3 c = mix(ambient, uBgColor, .3); //colour
         vec3 E = vec3(0.,0.,1.); //?
   
         // LOOP THROUGH LIGHT SOURCES
         
         for (int l = 0 ; l < nL ; l++) {
            // vec3 diff = uLc[l] * (diffuse * max(0., dot(norm, uLd[l])));
            // c+= diff;
            for(int n = 0; n < 8; n++){
               // COMPUTE NORMAL, INIT COLOR, APPROXIMATE VECTOR TO EYE
               vec3 N = normalize(P - S.xyz); //normal

            
               // ADD DIFFUSE AND SPECULAR ONLY IF NOT IN SHADOW
               

               float t = -1.;
               // float rayHalfspace(vec3 V, vec3 W, vec4 H)
               vec4 H = S * M;
               H /= sqrt(dot(H.xyz, H.xyz));
               t = max(t, rayHalfspace(P, uLd[l], H));
               N = H.xyz;
      
               // COMPUTE DIFFUSE AND SPECULAR FOR THIS LIGHT SOURCE
               //color += vec3(max(0., dot(Nt.xyz, vec3(.2))));
               if (t < 0.) {
                  vec3 R = 2. * dot(N, uLd[l]) * N - uLd[l];
                  vec3 reflectDir = reflect(-uLd[l], N);
                  c += uLc[l] * (
                     diffuse * max(0.,dot(N, uLd[l]))
                     + specular * pow(max(dot(E, reflectDir), 0.0), p)
                     );
               }
            }
         }
         //c *= 1. + .5 * noise(3.*N); // OPTIONAL SPOTTY TEXTURE
         return c;
      }
    
      void main() {
   
         // BACKGROUND COLOR IS THE DEFAULT COLOR
   
         vec3 color = uBgColor;
   
         // DEFINE RAY INTO SCENE FOR THIS PIXEL
   
         vec3 V = vec3(0.,0.,fl);
         vec3 W = normalize(vec3(vPos.xy, -fl));
   
         // LOOP THROUGH SPHERES
         float tMin = 10000.;
         for (int n = 0 ; n < nS ; n++) {
            float t = raySphere(V, W, uS[n]);
            if (t > 0. && t < tMin) {
               // IF THIS IS THE NEAREST SPHERE, DO PHONG SHADING
   
               vec3 P = V + t * W;
               color = shadeSphere(P, uS[n], uSm[n]);
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
                        rColor = shadeSphere(P + rt * R, uS[rn], uSm[rn]);
                     }
               }
               if (rtMin < 10000.)
                  color += .5 * rColor;
            }
         }
         
         // RAY TRACE TO THE POLY
         for (int n = 0; n < nP; n++){
            vec4 Nt = rayOctahedron(V, W, uIM[n]);
            vec3 P = V + Nt.w * W;
            
            if (Nt.w > 0. && Nt.w < tMin) {
            //shading poly
            color = shade(P, W, Nt.xyz, uSm[1]);

            //SEE WHETHER ANY OTHER SPHERE IS VISIBLE VIA REFLECTION
            for (int sn = 0; sn < 8; sn ++){
               vec3 N = normalize(P - uOctahedron[sn].xyz);
               vec3 R = 2. * dot(N, -W) * N + W;
               float rtMin = 10000.;
               vec3 rColor;

               //reflect spheres
               for (int rn = 0 ; rn < nS ; rn++) {
                  float rt = raySphere(P, R, uS[rn]);
                  if (rt > 0. && rt < rtMin) {
                     rtMin = rt;
                     rColor = shadeSphere(P + rt * R, uS[rn], uSm[rn]);
                  }
               }
               //reflect
               if (rtMin < 10000.)
                  color += .25 * rColor;
               }
               //color *= 1. + .5 * noise(3.* N); // OPTIONAL SPOTTY TEXTURE
            }
         }
         
   
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
       S.setUniform('Matrix4fv', 'uSm', false, S.material.flat());
    
       // SEND BACKGROUND COLOR TO GPU
    
       S.setUniform('3fv', 'uBgColor', [ red.value  /100,
                                         green.value/100,
                                         blue.value /100 ]);
    
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
         x,y,z,1]; //top -> down is right -> left

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
      for(let i = 0; i < S.nP; i++){
         let cM = matrixIdentity();
         cM = matrixMultiply(cM, matrixTranslate(.02*i, .04*i, 0));
         cM = matrixMultiply(cM, matrixRotz(time));
         // cM = matrixMultiply(cM, matrixRotx(time));
         cM = matrixMultiply(cM, matrixRoty(time));
         cM = matrixMultiply(cM, matrixTranslate(Math.cos(time * (.34 * i))/2, Math.sin(time * (.42 * i))/4, .2*i));
         
         cM = matrixMultiply(cM, matrixScale(.1, .15, .1));
         polyM.push(matrixInverse(cM));
      }
      S.setUniform('Matrix4fv', 'uIM', false, polyM.flat());
      
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
    
    
    