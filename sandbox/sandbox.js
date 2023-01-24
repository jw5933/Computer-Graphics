class Sandbox{
    constructor(x, y, cx = 16, cy = 32, cz = 16, ) {
        noise.seed(Math.random());
        let chunks = this.chunks = new Array(x);

        console.log("new sandbox");
        for(let i = 0; i < chunks.length; i++){
            chunks[i] = new Array(y);
            for(let j = 0; j < chunks[i].length; j++){
                chunks[i][j] = new Chunk(cz, cy, cz, {x:i,y:j});
            }
        }
        this.sx = x;
        this.sy = y;
        //                  dirt,   grass,  stone,  oak,    plank,  stone brick glass   orchid      leaf
        this.pickBlocks = {'1': 1,  '2': 2, '3': 3, '4': 4, '5': 5, '6': 8,     '7': 9, '8': 10,    '9': 11 }
        Chunk.prototype.totalChunks = {x: x, y: y};
    }

    createFlatBox = () => {
        console.log("create flat box", this.sx, this.sy);
        for (let x = 0; x < this.sx; x++){
            for(let y = 0; y < this.sy; y++){
                let height = Object.keys(BLOCK).length;
                this.chunks[x][y].createFlatChunk(height);
            }
        }
    }

    createWorld = () => {
        // console.log("create world");
        console.log("create world", this.sx, this.sy);
        for (let x = 0; x < this.sx; x++){
            for(let y = 0; y < this.sy; y++){
                this.chunks[x][y].createChunk();
            }
        }
    }
}