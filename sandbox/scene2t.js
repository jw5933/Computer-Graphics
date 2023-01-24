
rooms.scene2t = function() {

lib3D2();

description = `<b>Sandbox world</b>
               <p>
               3x3 chunks of 16x32x16 blocks. <br>If it doesn't load, please refresh or use incognito.
               <br>Click space if you lose highlighted block.
               <div id=selectInfo>selected block: dirt&nbsp;</div>
               <ul>
				   <li>WASD: move camera around in x and y.</li>
				   <li>Q/E: move camera in and out.</li>
				   <li>Drag canvas: look around.</li>
				   <li>arrow keys: pick block in x and z.</li>
				   <li>shift/space keys: pick block below/above in y.</li>
				   <li>I: insert selected block at highlighted block.</li>
				   <li>DEL: delete highlighted block.</li>
				   <li>P + KEY: place selected block in direction -> [<br>arrow keys: x & z, <br>shift: -y, <br>space: y] </li>
				   <li>selectable blocks (# keys): [<br>1: dirt,  <br>2: grass block, <br>3: stone, <br>4: wood, <br>5: wood plank, <br>6: stone brick, <br>7: glass, <br>8: flower, <br>9: leaf block ]</li>
               </ul>
				
               <br>
	       `;

code = {
'init':`





















	
	if (!S.sandbox){
		S.sSize = {x: 3, y: 3};
		S.cSize = {x: 16, y: 32, z: 16};
		S.sandbox = new Sandbox(S.sSize.x, S.sSize.y, S.cSize.x, S.cSize.y, S.cSize.z);
		
		S.sandbox.createWorld();
		// console.log(S.sandbox);
	}
	
	
	//movement
	S.pos = {x: -1, y: 0, z:-7};
	//key events
	S.pressedW = false;
	S.pressedA = false;
	S.pressedS = false;
	S.pressedD = false;
	S.pressedA = false;
	S.pressedLA= false;
	S.pressedRA = false;
	S.pressedUA = false;
	S.pressedDA = false;
	S.pressedDEL = false;
	S.pressedP = false;
	S.pressedI = false;
	S.pickBlock = 1;
	
	S.getBlock = (x, y, z) => {
		//returns chunk (sx,  sy), x, z, y, id
		let sx = Math.floor(x/S.cSize.x);
		let sy = Math.floor(z/S.cSize.z);
		x %=S.cSize.x;
		z %= S.cSize.z;
		let block = {sx, sy, x, z, y, id: S.sandbox.chunks[sx][sy].blocks[x][z][y]};
		// console.log(block);
		return block;
	}
	
	S.getTopBlock = (x, y, z) => {
		y = S.cSize.y-1;
		while(y >= 0 && S.getBlock(x, y, z).id == 0){
			y--;
		}
		return {x, y, z};
	}

	S.validateBlock = (x, y, z) => {
		// console.log("validating", x, y, z, x/S.cSize.x <= S.sSize.x && z/S.cSize.z <= S.sSize.y && y < S.cSize.y);
		return (x/S.cSize.x <= S.sSize.x && z/S.cSize.z <= S.sSize.y && y < S.cSize.y);
	}
	let validateBlock = S.validateBlock;
	
	S.visibleBlock = (x, y, z) => {
		let left = x-1 < 0 || S.getBlock(x-1, y, z).id == 0;
		let right = x+1 >= S.sSize.x * S.cSize.x || S.getBlock(x+1, y, z).id == 0;
		let top = y+1 >= S.cSize.y || S.getBlock(x, y+1, z).id == 0;
		let bot = y-1 < 0 || S.getBlock(x, y-1, z).id == 0;
		let back = z-1 < 0 || S.getBlock(x, y, z-1).id == 0;
		let front = z+1 >= S.sSize.y * S.cSize.z || S.getBlock(x, y, z+1).id == 0;
		return left||right||top||bot||back||front;
	}
	let visibleBlock = S.visibleBlock;
	
	S.getCloseBlock = (x, y, z, l = false, r = false, u = false, d = false, f= false, b = false, v = false) => {
		// console.log(x, y, z, l , r, u, d, f, b, v);
		let queue = [];
		queue.push({x: x, y: y, z: z});
		
		let visited = new Array(S.sSize.x * S.cSize.x);
		for(let i = 0; i < visited.length; i++){
			visited[i] = new Array(S.sSize.y * S.cSize.z);
			for(let j = 0; j < visited[0].length; j++){
				visited[i][j] = new Array(S.cSize.y);
			}
		}
		visited[x][z][y] = true;
		
		//returns closest block, always following down, up, right, back, left, ford
		while(queue.length > 0){
			let size = queue.length;
			for (let i = 0; i < size; i++) {
				let block = queue.shift();
				let x = block.x, y = block.y, z = block.z;
				//down
				if (d && block.y - 1 >= 0 && validateBlock(x, y-1, z) && !visited[x][z][y-1]) {
					block = {x: x, y: y-1, z: z};
					if (S.getBlock(x, y-1, z).id != 0 && (!v || visibleBlock(x, y-1, z))) 
						return block;
					queue.push(block);
					visited[x][z][y-1] = true;
				}
				//up
				if (u && validateBlock(x, y+1, z) && !visited[x][z][y+1]) {
					block = {x: x, y: y+1, z: z};
					if (S.getBlock(x, y+1, z).id != 0 && (!v ||visibleBlock(x, y+1, z))) 
						return block;
					queue.push(block);
					visited[x][z][y+1] = true;
				}
				//right
				if (r && validateBlock(x+1, y, z) && !visited[x+1][z][y]) {
					block = {x: x+1, y: y, z: z};
					console.log(block);
					if (S.getBlock(x+1, y, z).id != 0 && (!v ||visibleBlock(x+1, y, z))) 
						return block;
					queue.push(block);
					visited[x+1][z][y] = true;
				}
				//back
				if (b && block.z - 1 >= 0 && validateBlock(x, y, z-1) && !visited[x][z-1][y]) {
					block = {x: x, y: y, z: z-1};
					if (S.getBlock(x, y, z-1).id != 0 && (!v ||visibleBlock(x, y, z-1))) 
						return block;
					queue.push(block);
					visited[x][z-1][y] = true;
				}
				//left
				if (l && x-1 >= 0 && validateBlock(x-1, y, z) && !visited[x-1][z][y]) {
					block = {x: x-1, y: y, z: z};
					if (S.getBlock(x-1, y, z).id != 0 && (!v ||visibleBlock(x-1, y, z))) 
						return block;
					queue.push(block);
					visited[x-1][z][y] = true;
				}
				//ford
				if (f && validateBlock(x, y, z+1) && !visited[x][z+1][y]) {
					block = {x: x, y: y, z: z+1};
					if (S.getBlock(x, y, z+1).id != 0 && (!v ||visibleBlock(x, y, z+1))) 
						return block;
					queue.push(block);
					visited[x][z+1][y] = true;
				}
			}
		}
		console.log("no block found");
		return null; //no blocks found
	}
	
	S.setBlock = (block) => {
		if (block){
			if (!S.visibleBlock(block.x, block.y, block.z)) 
				block = S.getTopBlock(block.x, block.y, block.z);
			
			S.currentBlock = block;
		}
	}
	
	S.currentBlock = S.getCloseBlock(0, 0, 0, true, true, true, true, true, true, true);
	S.currentBlock = S.getTopBlock(S.currentBlock.x, S.currentBlock.y, S.currentBlock.z);
	if(S.currentBlock!= null) console.log(S.currentBlock);
	
	S.rot = {x: 0, y: 0};
	S.initialRot = {x:0, y:0};
	
	S.calculateAngle = (x, y) => {
		let i = S.initialRot;
		let r = S.rot;
		r.x = (r.x - (i.x - x) + 360) % 360;
		r.y = clamp(-90, 90, r.y + (i.y - y));
		// console.log(r.x, r.y);
	}
	
	let clamp = (min, max, val) => {
		let r = Math.min(max, Math.max(val, min));
		return r;
	}
	
	S.material = [
		[.1,.1,.1,0,     .9,.9,.9,0,  0,0,0,5,    0,0,0,0], // PAPER
		[.1,.1,.1,0,     .2,.2,.2,0,  1,1,1,5,    0,0,0,0], // SILVER
		[.25,0,0,0,      .5,0,0,0,    2,2,2,20,   0,0,0,0], // RED PLASTIC
		[.15,.05,.025,0, .3,.1,.05,0, .6,.2,.1,3, 0,0,0,0], // COPPER
		[.25,.15,.025,0, .5,.3,.05,0, 1,.6,.1,6,  0,0,0,0], // GOLD
		[.05,.05,.05,0,  .1,.1,.1,0,  1,1,1,5,    0,0,0,0], // LEAD
		];
	S.nM = S.material.length;

   // TWO TRIANGLES THAT ARE NOT A TRIANGLE STRIP

   S.twoTrianglesMesh = [
      0,0,0,   0,0,1,  0,0,
      1,0,0,   0,0,1,  1,0,
      0,1,0,   0,0,1,  0,1,

      0,0,0,   0,1,0,  0,0,
      0,0,1,   0,1,0,  1,0,
      1,0,0,   0,1,0,  0,1,
   ];

   S.twoTrianglesMesh.isTriangles = true;

   // A SQUARE IS A TRIANGLE MESH WITH JUST TWO TRIANGLES

   S.squareMesh = [ -1, 1, 0,  0,0,1,  0,1,
   					-1,-1, 0,  0,0,1,  0,0,
                     1, 1, 0,  0,0,1,  1,1,
                     1,-1, 0,  0,0,1,  1,0];
                     
   S.squareMesh2 = [ 0, 0, 0, 0, 0, 1, 0.5, 0.5,
   				     1, 1, 0,  0,0,1,  1,1,
   				     -1, 1, 0,  0,0,1,  0,1,
   				     -1, -1, 0,  0,0,1,  0,0,
   				     1, -1, 0,  0,0,1,  1,0,
   				      1, 1, 0,  0,0,1,  1,1];

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

   let uvMesh = (f,nu,nv,data) => {
      let mesh = [];
      for (let iv = 0 ; iv < nv ; iv++) {
         let v = iv / nv;
	 let strip = [];
         for (let iu = 0 ; iu <= nu ; iu++) {
	    let u = iu / nu;
	    strip = strip.concat(f(u,v,data));
	    strip = strip.concat(f(u,v+1/nv,data));
	 }
	 mesh = glueMeshes(mesh, strip);
      }
      return mesh;
   }

   S.uvMesh = uvMesh;

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
   let face = transformMesh(S.squareMesh,  matrixRotz( Math.PI/2));
   // face = S.squareMesh;
   face = transformMesh(face, 			   matrixTranslate([0,0,1]));
   let face0 = transformMesh(face, 		   matrixRotz( Math.PI/2));
   let face1 = transformMesh(face0,        matrixRoty( Math.PI/2));
   let face2 = transformMesh(face1,        matrixRoty( Math.PI/2));
   let face3 = transformMesh(face2,        matrixRoty(Math.PI/2));
   let face4 = transformMesh(face0,        matrixRotx(-Math.PI/2));
   let face5 = transformMesh(face0,        matrixRotx( Math.PI/2));
   S.cubeMesh = glueMeshes(face0,
                glueMeshes(face1,
                glueMeshes(face2,
                glueMeshes(face3,
                glueMeshes(face4,
		           face5)))));
		           
   face = transformMesh(S.squareMesh2,  matrixRotz( Math.PI/2))
   // face = transformMesh(face, 			   matrixTranslate([0,0,0.5]));
   face0 = transformMesh(face, 		   matrixRotz( Math.PI/2));
   face1 = transformMesh(face0,        matrixRoty( Math.PI/2));
   S.plantMesh = [face0, [0, 0, 0, 0, 0, 1, 0.5, 0.5], face1].flat();
   
   // [face0, face1, face2, face3, [0, 0, 0, 0, 0, 1, 0.5, 0.5]].flat();
   
   // glueMeshes(face0, 
   // 				 glueMeshes(face1, 
   // 				 glueMeshes(face2, 
   // 				 face3)));

   S.applyTexture = (mesh, texture, start, count) => {
      let gl = S.gl;
      if (texture) {
         if (! S.textures[texture]) { // NEED TO LOAD THE TEXTURE FROM THE SERVER.
            let image = new Image();
            image.onload = function() {
               S.textures[this.texture] = gl.createTexture();
               gl.bindTexture     (gl.TEXTURE_2D, S.textures[this.texture]);
               gl.texImage2D      (gl.TEXTURE_2D, 0, gl.RGBA, gl.RGBA, gl.UNSIGNED_BYTE, this);
               gl.texParameteri   (gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
               gl.texParameteri   (gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR_MIPMAP_NEAREST);
               gl.generateMipmap  (gl.TEXTURE_2D);
            }
            image.texture = texture;
            image.src = texture;
         }
         else {                          // TEXTURE HAS LOADED. OK TO RENDER.
            gl.activeTexture(gl.TEXTURE0);
            gl.bindTexture(gl.TEXTURE_2D, S.textures[texture]);
         }
      }

      gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(mesh), gl.STATIC_DRAW);
      gl.drawArrays(mesh.isTriangles ? S.gl.TRIANGLES
                                     : S.gl.TRIANGLE_STRIP, start, count);
   }

   // DRAW A SINGLE MESH. WE STILL NEED TO ADD MATERIAL PROPERTIES!

   S.textures = {};
   S.drawMesh = (mesh, matrix, materialIndex, textureData, highlight) => {
		let gl = S.gl;
		S.setUniform('Matrix4fv', 'uMatrix', false, matrix);
		S.setUniform('Matrix4fv', 'uInvMatrix', false, matrixInverse(matrix));
		S.setUniform('Matrix4fv', 'uMaterial', false, S.material[materialIndex]);

		S.setUniform('1i', 'uSampler', 0);
		S.setUniform('1f', 'uTexture', textureData ? 1 : 0);
		S.setUniform('3fv', 'uHighlight', highlight ? [1,0,0] : [1,1,1]);
		
		let start = 0;
		if (textureData == null){
			gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(mesh), gl.STATIC_DRAW);
			gl.drawArrays(mesh.isTriangles ? S.gl.TRIANGLES
										 : S.gl.TRIANGLE_STRIP, 0, mesh.length / S.VERTEX_SIZE);
			return;
	 	}
		
		for(let i = 0; i < textureData.length; i++){
			let count = textureData[i].count;
			S.applyTexture(mesh, textureData[i].texture, start, count);
			start += count;
		}
   }
   
   //matrix, id -> put type, texture data into block data
	S.createBlock = (m, id, highlight) => {
		//get block data? by id
		let block = BLOCKFUNC.getBlock(id);
		
		//block function
		mesh = S[block.type];
		S.drawMesh(mesh, m.get(), 0, block.textureData(), highlight);
	}
	
	S.drawScene = () => {
		S.setUniform('3fv', 'uCamPos', [S.pos.x, S.pos.y, S.pos.z]);
	   // SET THE PROJECTION MATRIX BASED ON CAMERA FOCAL LENGTH
	
	   let fl = 5.0;
	   let value = Math.max(0, 1 + S.pos.z/10);
	   S.setUniform('Matrix4fv', 'uProject', false,
		  [value,0,0,0, 0,value,0,0, 0,0,1,-1/fl, 0,0,0,value]);
	
	   // SPECIFY SCENE LIGHTING
	   S.nL = 2;
	   S.setUniform('3fv', 'uLd', [ .57,.57,.57, -.57,-.57,-.57 ]);
	   S.setUniform('3fv', 'uLc', [ 1,1,1, .5,.3,.1 ]);
	   S.setUniform('3fv', 'uBgColor', [ .89,.81,.75 ]);
	   
	   
	   // RENDER THE SCENE
	   let textureData = [];
	   let blockScale = 0.1;
	   let blockSize = blockScale * 2;
	   /////
	   let m = new Matrix();
	   let cBlock = S.getBlock(S.currentBlock.x, S.currentBlock.y, S.currentBlock.z);
	   // console.log(cBlock);
		m.identity();
		m.save();
			m.translate(-0.5, -1, 0);
			m.rotx(Math.PI/4 + S.rot.y/180 * Math.PI);
			m.roty(Math.PI/2 + S.rot.x/180 * Math.PI);
	
			m.save();
			let sandbox = S.sandbox;
			for(let sx = sandbox.sx-1; sx >= 0; sx--){ //x
			m.save();
				for (let sy = sandbox.sy-1; sy >=0; sy--){
				m.save();
					let chunk = sandbox.chunks[sx][sy];
					m.translate(sx * chunk.sx * blockSize, 0, sy * chunk.sz * blockSize);
	
					for (let x = chunk.sx-1; x >=0; x--){
					m.save();
						m.translate(x * blockSize, 0, 0);
						let l = x - 1;
						let r = x + 1;
						
						let lChunk = chunk;
						if (l < 0){
							l = chunk.sx-1;
							lChunk = sx-1 >= 0 ? sandbox.chunks[sx-1][sy] : {blocks: chunk.emptyBlock(chunk.sx, chunk.sy, chunk.sz)};
						}
	
						let rChunk = chunk;
						if (r >= chunk.sx){
							r = 0;
							rChunk = sx+1 < sandbox.sx ? sandbox.chunks[sx+1][sy] : {blocks: chunk.emptyBlock(chunk.sx, chunk.sy, chunk.sz)};
						}
						
						for (let z = chunk.sz-1; z >=0; z--){
						m.save();
							m.translate(0, 0, z * blockSize);
							let f = z + 1;
							let b = z - 1;
							let bChunk = chunk;
							if (b < 0){
								b = chunk.sx-1;
								bChunk = sy-1 >= 0 ? sandbox.chunks[sx][sy-1] : {blocks: chunk.emptyBlock(chunk.sx, chunk.sy, chunk.sz)};
							}
		
							let fChunk = chunk;
							if (f >= chunk.sx){
								f = 0;
								fChunk = sy+1 < sandbox.sy ? sandbox.chunks[sx][sy+1] : {blocks: chunk.emptyBlock(chunk.sx, chunk.sy, chunk.sz)};
							}
							
							for (let y = 0; y < S.cSize.y; y++){
							m.save();
								let blockId = chunk.blocks[x][z][y];
								// console.log(blockId, top);
								if (blockId == 0) {
									m.restore();
									continue; //ok many subsequent blocks air or not air
								}
								
								let left = lChunk.blocks[l][z][y];
								let right = rChunk.blocks[r][z][y];
								let ford = fChunk.blocks[x][f][y];
								let back = bChunk.blocks[x][b][y];
								let top = y + 1 >= S.cSize.y ? 0: chunk.blocks[x][z][y+1];
								
								//only render those that are seen
								if (top != 0 && left != 0 && right != 0 && ford != 0 && back != 0){ 
									top = blockId;
									m.restore();
									continue;
								}
								let highlight = cBlock.sx == sx && cBlock.sy == sy 
									&& cBlock.x == x && cBlock.y == y && cBlock.z == z;
									
								m.translate(0, y * blockSize, 0);
								m.scale(blockScale);
								
								S.createBlock(m, blockId, highlight);
								top = blockId;
							m.restore();
							}
						m.restore();
						}
					m.restore();
					}
				m.restore();
				}
			m.restore();
			}
		  m.restore();
		  
	   m.restore();
   }
`,
vertex: `
S.setVertexShader(\`
   uniform vec3 uCamPos;
   attribute vec3 aPos, aNor;
   attribute vec2 aUV;
   uniform   mat4 uMatrix, uInvMatrix, uProject;
   varying   vec3 vPos, vNor;
   varying   vec2 vUV;

   void main() {
      vPos = (uProject * uMatrix * vec4(aPos, 1.)).xyz + uCamPos;
      vNor = (vec4(aNor, 0.) * uInvMatrix).xyz;
      vUV = aUV;
      vec3 pos = vec3(vPos.xy,-.01 * vPos.z);
      gl_Position = vec4(pos, 1.);
   }
\`)
`,
fragment: `
S.setFragmentShader(\`
   const int nL = \` + S.nL + \`;
   const int nM = \` + S.nM + \`;
   uniform vec3 uBgColor;
   uniform vec3 uLd[nL];
   uniform vec3 uLc[nL];
   uniform mat4 uMaterial;
   uniform vec3 uHighlight;

   uniform sampler2D uSampler; // INDEX OF THE TEXTURE TO BE SAMPLED
   uniform float uTexture;     // ARE WE RENDERING TEXTURE FOR THIS OBJECT?

   varying vec3 vPos, vNor;
   varying vec2 vUV;

   void main() {
      vec3 N = normalize(vNor);
      vec3  ambient  = uMaterial[0].rgb;
      vec3  diffuse  = uMaterial[1].rgb;
      vec3  specular = uMaterial[2].rgb;
      float p        = uMaterial[2].a;
      vec3 c = mix(ambient, uBgColor, .3);
      for (int l = 0 ; l < nL ; l++) {
         vec3 R = 2. * dot(N, uLd[l]) * N - uLd[l];
         c += uLc[l] * (diffuse * max(0.,dot(N, uLd[l]))
                      + specular * pow(max(0., R.z), p));
      }

      vec4 texture = texture2D(uSampler, vUV);
      float t = vNor.z;
      // c = mix(vec3(1.5,1.5,1.5), vec3(1.,-.5,-2.), t);
      c *= mix(vec3(1.), texture.rgb, texture.a * uTexture);
      float alpha = mix(texture.a, 1., 1. - uHighlight.g);
	  gl_FragColor = (alpha) * vec4(c, 1.) * vec4(uHighlight, 1.);
      

	  // gl_FragColor = texture;
      // gl_FragColor = vec4(c, 1.);
   }
   \`);
// S.setFragmentShader2(\`
//    const int nL = \` + S.nL + \`;
//    const int nM = \` + S.nM + \`;
//    uniform vec3 uBgColor;
//    uniform vec3 uLd[nL];
//    uniform vec3 uLc[nL];
//    uniform mat4 uMaterial;
//    uniform float uPick;
//
//    uniform sampler2D uSampler; // INDEX OF THE TEXTURE TO BE SAMPLED
//    uniform float uTexture;     // ARE WE RENDERING TEXTURE FOR THIS OBJECT?
//
//    varying vec3 vPos, vNor;
//    varying vec2 vUV;
//
//    void main() {
//       vec3 N = normalize(vNor);
//       vec3  ambient  = uMaterial[0].rgb;
//       vec3  diffuse  = uMaterial[1].rgb;
//       vec3  specular = uMaterial[2].rgb;
//       float p        = uMaterial[2].a;

//       vec4 texture = texture2D(uSampler, vUV);
//       gl_FragColor = texture;
// }
// \`);
`,
render: `
	//movement
	if (S.pressedW)
		S.pos.y -= 0.1;
	if (S.pressedA)
		S.pos.x += 0.1;
	if(S.pressedS)
		S.pos.y += 0.1;
	if(S.pressedD)
		S.pos.x -= 0.1;
	if (S.pressedQ)
		S.pos.z += 0.3;
	if (S.pressedE)
		S.pos.z -= 0.3;
	
	S.drawScene();
	
`,
events: `
	onKeyPress   = key => {
		S.pressedW = key == 87 || S.pressedW;
		S.pressedA = key == 65 || S.pressedA;
		S.pressedS = key == 83 || S.pressedS;
		S.pressedD = key == 68 || S.pressedD;
		S.pressedQ = key == 81 || S.pressedQ;
		S.pressedE = key == 69 || S.pressedE;
		S.pressedShift = key == 16 || S.pressedShift;
		S.pressedSpace = key == 32 || S.pressedSpace;
		S.pressedLA = key == 37 || S.pressedLA;
		S.pressedRA = key == 39 || S.pressedRA;
		S.pressedUA = key == 38 || S.pressedUA;
		S.pressedDA = key == 40 || S.pressedDA;
		// S.pressedDEL = key == 8 || S.pressedDEL;
		S.pressedP = key == 80 || S.pressedP;
		// S.pressedI = key == 73 || S.pressedI;
		let block = S.currentBlock;
		
		//navigation
		let b = null;
		switch(key){
		 case 40: //down arrow
		 	console.log("down arrow");
			if (block.x - 1 >= 0){
				if (S.pressedP){
					b = S.getBlock(block.x-1, block.y, block.z);
					S.sandbox.chunks[b.sx][b.sy].blocks[b.x][b.z][b.y] = S.pickBlock;
					return;
				}
				let newBlock = S.getTopBlock(block.x-1, block.y, block.z);
				S.setBlock(newBlock);
			}
		 break;
		 case 38: //up arrow
		 	console.log("up arrow");
			if (block.x + 1 < S.sSize.x * S.cSize.x){
				if (S.pressedP){
					b = S.getBlock(block.x+1, block.y, block.z);
					S.sandbox.chunks[b.sx][b.sy].blocks[b.x][b.z][b.y] = S.pickBlock;
					return;
				}
				let newBlock = S.getTopBlock(block.x+1, block.y, block.z);
				S.setBlock(newBlock);
			}
		 break;
		 case 39: //right arrow
		 	console.log("right arrow");
			if (block.z + 1 < S.sSize.y * S.cSize.z){
				if (S.pressedP){
					b = S.getBlock(block.x, block.y, block.z+1);
					S.sandbox.chunks[b.sx][b.sy].blocks[b.x][b.z][b.y] = S.pickBlock;
					return;
				}
				let newBlock = S.getTopBlock(block.x, block.y, block.z+1);
				S.setBlock(newBlock);
			}
		 break;
		 case 37: //left arrow
		 	console.log("left arrow");
			if (block.z - 1 >= 0){ 
				if (S.pressedP){
					b = S.getBlock(block.x, block.y, block.z-1);
					S.sandbox.chunks[b.sx][b.sy].blocks[b.x][b.z][b.y] = S.pickBlock;
					return;
				}
				let newBlock = S.getTopBlock(block.x, block.y, block.z-1);
				S.setBlock(newBlock);
			}
		 break;
		 case 16: //shift
		 	if (block.y - 1 >= 0){ 
		 		if (S.pressedP){
					b = S.getBlock(block.x, block.y-1, block.z);
					S.sandbox.chunks[b.sx][b.sy].blocks[b.x][b.z][b.y] = S.pickBlock;
					return;
				}
				let newBlock = S.getCloseBlock(block.x, block.y, block.z, undefined, undefined, undefined, true, undefined, undefined, true);
				S.setBlock(newBlock);
			}
		 break;
		 case 32: //space
		 	if (block.y + 1 < S.cSize.y){ 
		 		if (S.pressedP){
					b = S.getBlock(block.x, block.y+1, block.z);
					S.sandbox.chunks[b.sx][b.sy].blocks[b.x][b.z][b.y] = S.pickBlock;
					return;
				}
				let newBlock = S.getCloseBlock(block.x, block.y, block.z, undefined, undefined, true, undefined, undefined, undefined, true);
				S.setBlock(newBlock);
			}
		 break;
		 case 8: //delete
		 	b = S.getBlock(block.x, block.y, block.z);
			S.sandbox.chunks[b.sx][b.sy].blocks[b.x][b.z][b.y] = 0;
			S.setBlock(S.getCloseBlock(block.x, block.y, block.z, true, true, true, true, true, true, true));
		 break;
		 case 73: //insert (replace)
		 	b = S.getBlock(block.x, block.y, block.z);
			S.sandbox.chunks[b.sx][b.sy].blocks[b.x][b.z][b.y] = S.pickBlock;
		 break;
		 	
		 break;
		 case 49: //1
		 case 50: //2
		 case 51: //3
		 case 52: //4
		 case 53: //5
		 case 54: //6
		 case 55: //7
		 case 56: //8
		 case 57: //9
		 	let v = (key - 48);
		 	console.log(v, S.sandbox.pickBlocks[v]);
		 	S.pickBlock = S.sandbox.pickBlocks[v];
		 	selectInfo.innerHTML = 'selected block: ' + BLOCKFUNC.getBlockName(S.sandbox.pickBlocks[v]);
		 break;
		}
	}
	
	
  	onKeyRelease = key => {
  		S.pressedW = !(key == 87) && S.pressedW;
  		S.pressedA = !(key == 65) && S.pressedA;
  		S.pressedS = !(key == 83) && S.pressedS;
  		S.pressedD = !(key == 68) && S.pressedD;
  		S.pressedQ = !(key == 81) && S.pressedQ;
  		S.pressedE = !(key == 69) && S.pressedE;
  		S.pressedShift = !(key == 16) && S.pressedShift;
		S.pressedSpace = !(key == 32) && S.pressedSpace;
		S.pressedLA= !(key == 37) && S.pressedLA;
		S.pressedRA = !(key == 39) && S.pressedRA;
		S.pressedUA = !(key == 38) && S.pressedUA;
		S.pressedDA = !(key == 40) && S.pressedDA;
		// S.pressedDEL = !(key == 8) && S.pressedDEL;
		S.pressedP = !(key == 80) && S.pressedP;
		// S.pressedI = !(key == 73) && S.pressedI;
  	}
  	
  	onPress = (x, y) => {
  		S.initialRot.x = x;
  		S.initialRot.y = y;
  	}
  	
  	onDrag = (x,y) => {
     S.calculateAngle(x, y);
  }
   
`
};

}

