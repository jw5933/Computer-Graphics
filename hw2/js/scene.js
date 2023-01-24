
rooms.scene = function() {

    lib3D();
    
    description = `Raytracing a gold fish <br>in a fragment shader
    <small>
        <p>
        <b>Background color</b>
        <br> <input type=range id=red   value= 5> red
        <br> <input type=range id=green value=10> green
        <br> <input type=range id=blue  value=50> blue
        </p>

        <p>
        <b>Speed</b>
        <br> <input type = range id=timeScale value =1 min = 0.5 max = 5 step = 0.5> faster?
        </p>
    </small>
    `;
    
    code = {
    'init':`
    
       // DEFINE NUMBER OF SPHERES AND NUMBER OF LIGHTS
    
       S.nS = 100;
       S.nP = 5;
       S.nL = 2;
    
       // DEFINE MATERIALS TO BE RENDERED VIA PHONG REFLECTANCE MODEL
    
       let materials = [
          [.15,.05,.025,0, .3,.1,.05,0, .6,.2,.1,3, 0,0,0,0], // COPPER
          [.25,.15,.025,0, .5,.3,.05,0, 1,.6,.1,6,  0,0,0,0], // GOLD
          [.25,0,0,0,      0,0,0,0,    2,2,2,20,   0,0,0,0], // PLASTIC
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
       uniform vec4 uS;
       uniform mat4 uSm[nS];
       uniform mat4 uIM[5];
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
            //no shadows
   
            // COMPUTE DIFFUSE AND SPECULAR FOR THIS LIGHT SOURCE
   
            if (t < 0.) {
               vec3 R = 2. * dot(N, uLd[l]) * N - uLd[l];
               c += uLc[l] * (diffuse * max(0.,dot(N, uLd[l]))
                           + specular * pow(max(0., dot(R, E)), p));
            }
         }
          return c;
       }

      // SHADE A POLY
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
            for(int n = 1; n < 4; n ++){
                t = max(t, rayPyramid(P, uLd[l], uIM[n]).w);
            }

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

        // RENDER IN ORDER
        float tMin = 10000.;
        
        // RAY TRACE TO THE POLY
        vec4 Nt = vec4(0.0, 0.0, 0.0, 0.0);
        for (int n = 0; n < 4; n ++){
            if (n == 0){
                Nt = rayOctahedron(V, W, uIM[n]);
            }
            else{
                Nt = rayPyramid(V, W, uIM[n]);
            }

            vec3 P = V + Nt.w * W;
            if (Nt.w > 0. && Nt.w < tMin) {
                //shading poly
                color = shade(P, W, Nt.xyz, uSm[1]);
            }

            //no reflections to not break the illusion of fish..
        }

        float t = raySphere(V, W, uS);
        if (t > 0. && t < tMin) {
            // IF THIS IS THE NEAREST SPHERE, DO PHONG SHADING
            vec3 P = V + t * W;
            color = shadeSphere(P, uS, uSm[2]);
            tMin = t;
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

      const TAU = 6.2831853071;

      // ====================== DATA =========================//
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

      // SEND SPHERES DATA TO GPU
      let sData = [];
      for (let n = 0 ; n < S.nS ; n++)
         sData.push(S.sPos[n], radius);
        
      S.setUniform('Matrix4fv', 'uSm', false, S.material.flat());
   
      // SEND BACKGROUND COLOR TO GPU
      S.setUniform('3fv', 'uBgColor', [ red.value  /100,
                                       green.value/100,
                                       blue.value /100 ]);
   
   
      // Creating multiple polys using translations

      let newTime = time * timeScale.value;

      let polyM = [];
      let iM = matrixIdentity();

      //final movement (6)
      let finalM = iM;
      finalM = matrixMultiply(finalM, matrixRotz(Math.sin(newTime * TAU * 0.5) * 0.3));

      //body (0)
      let body = iM;
      let bodySize = 0.45;
      body = matrixMultiply(body, finalM);
      body = matrixMultiply(body, matrixTranslate(0, 0, 0.2));
      body = matrixMultiply(body, matrixScale(bodySize, bodySize, bodySize));
      //polyM.push(matrixInverse(body));


      //tail (1)
      let tail = iM;
      let tailSize = .2;
      
      tail = matrixMultiply(tail, finalM);
      tail = matrixMultiply(tail, matrixTranslate(-0.6, 0, 0));

      tail = matrixMultiply(tail, matrixTranslate(0.6, 0, 0));
      tail = matrixMultiply(tail, matrixRoty(Math.sin(newTime * TAU)*0.45));
      tail = matrixMultiply(tail, matrixTranslate(-0.6, 0, 0));
      
      tail = matrixMultiply(tail, matrixRotx(TAU/4));
      tail = matrixMultiply(tail, matrixRotz(-TAU/4));
      tail = matrixMultiply(tail, matrixScale(tailSize*.15, tailSize*1.2, tailSize));

      //polyM.push(matrixInverse(tail));
      

      //fin 1
      let fin = iM;
      let finSize = .1;
      
      fin = matrixMultiply(fin, finalM);
      fin = matrixMultiply(fin, matrixTranslate(0.1, -0.1, 0.5));

      fin = matrixMultiply(fin, matrixTranslate(0.1, -0.1, 0));
      fin = matrixMultiply(fin, matrixRoty(Math.sin(newTime * TAU)*0.45));
      fin = matrixMultiply(fin, matrixTranslate(-0.1, 0.1, 0));
      
      fin = matrixMultiply(fin, matrixTranslate(0.1, -0.1, 0));
      fin = matrixMultiply(fin, matrixRoty(0.5));
      fin = matrixMultiply(fin, matrixTranslate(-0.1, 0.1, 0));

      fin = matrixMultiply(fin, matrixRotz(0.2));
      fin = matrixMultiply(fin, matrixRotx(TAU/4));
      fin = matrixMultiply(fin, matrixRotz(-TAU/4));

      fin = matrixMultiply(fin, matrixScale(finSize* 0.1, finSize*1.2, finSize));
      //polyM.push(matrixInverse(fin));

      //fin 2
      let fin2 = iM;

      fin2 = matrixMultiply(fin2, finalM);
      
      fin2 = matrixMultiply(fin2, matrixTranslate(0.1, -0.23, 0.5));

      fin2 = matrixMultiply(fin2, matrixTranslate(0.1, -0.1, 0));
      fin2 = matrixMultiply(fin2, matrixRoty(Math.sin(newTime * TAU)*0.45));
      fin2 = matrixMultiply(fin2, matrixTranslate(-0.1, 0.1, 0));

      fin2 = matrixMultiply(fin2, matrixRotz(1));
      fin2 = matrixMultiply(fin2, matrixRotx(TAU/4));
      fin2 = matrixMultiply(fin2, matrixRotz(-TAU/4));

      
      fin2 = matrixMultiply(fin2, matrixScale(finSize* 0.1, finSize*1.5, finSize*0.5));

      //eyes (5)
      let eye = iM;
      let eyeSize = 0.01;
      eye = matrixMultiply(eye, matrixScale(eyeSize, eyeSize, eyeSize));

      let eyeData = [0.28, Math.sin(newTime * TAU * 0.5)* 0.1, 0, 0.03];
      S.setUniform('4fv', 'uS', eyeData);

      //tell gpu
      polyM.push(matrixInverse(body));
      polyM.push(matrixInverse(tail));
      polyM.push(matrixInverse(fin2));
      polyM.push(matrixInverse(fin));
      polyM.push(matrixInverse(eye));

      
      S.setUniform('Matrix4fv', 'uIM', false, polyM.flat());

      //shapes
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
    
    
    