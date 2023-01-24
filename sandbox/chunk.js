class Chunk{
    totalChunks = [0,0];
    constructor(x, y, z, chunk) {
        console.log("new chunk")
        let blocks = this.blocks = new Array(x);
        for(let i = 0; i < blocks.length; i++){
            blocks[i] = new Array(z);
            for(let j = 0; j < blocks[i].length; j++){
                blocks[i][j] = new Array(y);
            }
        }
        this.sx = x;
        this.sy = y;
        this.sz = z;
        this.chunk = chunk;
    }

    createFlatChunk = (height) => {
        console.log(height, this.sx, this.sz, this.sy);
        for ( let x = 0; x < this.sx; x++ ) {
            for (let z = 0; z < this.sz; z++) {
                for (let y = 0; y < this.sy; y++) {
                    let size = Object.keys(BLOCK).length-1;
                    if (y < height) this.blocks[x][z][y] = y%size+1;
                }
            }
        }
        // this.blocks[0][0][height-1] = 0;
        // this.blocks[this.sx-1][this.sz-1][height-2] = 0;
    }

    createChunk = (amp = 20) => {
        let startx = this.sx * this.chunk['x'];
        let startz = this.sz * this.chunk['y'];
        // let totalSize = {x: Chunk.prototype.totalChunks['x'] * this.sx, y: Chunk.prototype.totalChunks['y'] * this.sz};
        // console.log(startx,startz, totalSize);
        for ( let x = 0; x < this.sx; x++ ) {
            for (let z = 0; z < this.sz; z++) {
                let height = Math.floor(Math.max(1, this.sy/2 + noise.fractal((startx+x)/amp, (startz+z)/amp)*(this.sy/4)));
                // console.log(height);
                for (let y = height; y < this.sy; y++) {
                    this.blocks[x][z][y] = 0;
                }
                for (let y = 0; y < height; y++) {
                    this.blocks[x][z][y] = 3;
                }
                for (let y = Math.max(0, height-3); y < height-1; y++) {
                    this.blocks[x][z][y] = 1;
                }
                this.blocks[x][z][Math.max(0, height-1)] = 2;
            }
        }
    }

    emptyBlock = (x, y, z) => {
        // console.log("creating empty blocks");

        let blocks = new Array(x);
        for(let i = 0; i < blocks.length; i++){
            blocks[i] = new Array(z);
            for(let j = 0; j < blocks[i].length; j++){
                blocks[i][j] = new Array(y);
                for(let k = 0; k < y; k++){
                    blocks[i][j][k] = 0;
                }
            }
        }
        // console.log(blocks);
        return blocks;
    }
}