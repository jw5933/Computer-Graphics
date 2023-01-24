//all textures via resource pack for minecraft: https://faithful.team
BLOCK = {};
BLOCKFUNC = {};

// air
BLOCK.air = {
    id: 0,
    type: "none",
};

BLOCK.dirt = {
    name: 'dirt',
    id: 1,
    type: "cubeMesh",
    textureData: () =>
    {
        let textureData = [];
        textureData.push({texture: 'textures/dirt.png', count: S.VERTEX_SIZE * 5 + 2}); //??
        return textureData;
    }
};

BLOCK.grass = {
    name: 'grass',
    id: 2,
    type: "cubeMesh",
    textureData: () =>
    {
        let textureData = [];
        textureData.push({texture: 'textures/grass_block_side.png', count: S.VERTEX_SIZE * 3}); //??
        textureData.push({texture: 'textures/grass_block_top.png', count: S.VERTEX_SIZE/2+1});
        textureData.push({texture: 'textures/dirt.png', count: S.VERTEX_SIZE/2+1});
        return textureData;
    }
};

BLOCK.stone = {
    name: 'stone',
    id: 3,
    type: "cubeMesh",
    textureData: () =>
    {
        let textureData = [];
        textureData.push({texture: 'textures/stone.png', count: S.VERTEX_SIZE * 5});
        return textureData;
    }
};

BLOCK.oakWood = {
    name: 'wood',
    id: 4,
    type: "cubeMesh",
    textureData: () =>
    {
        let textureData = [];
        textureData.push({texture: 'textures/oak_log.png', count: S.VERTEX_SIZE * 3}); //??
        textureData.push({texture: 'textures/oak_log_top.png', count: S.VERTEX_SIZE +2});
        return textureData;
    }
};

BLOCK.oakPlank = {
    name: 'wood plank',
    id: 5,
    type: "cubeMesh",
    textureData: () =>
    {
        let textureData = [];
        textureData.push({texture: 'textures/oak_planks.png', count: S.VERTEX_SIZE * 5});
        return textureData;
    }
};

BLOCK.spruceWood = {
    name: 'spruce wood',
    id: 6,
    type: "cubeMesh",
    textureData: () =>
    {
        let textureData = [];
        textureData.push({texture: 'textures/spruce_log.png', count: S.VERTEX_SIZE * 3}); //??
        textureData.push({texture: 'textures/spruce_log_top.png', count: S.VERTEX_SIZE +2});
        return textureData;
    }
};

BLOCK.sprucePlank = {
    name: 'spruce plank',
    id: 7,
    type: "cubeMesh",
    textureData: () =>
    {
        let textureData = [];
        textureData.push({texture: 'textures/spruce_planks.png', count: S.VERTEX_SIZE * 5});
        return textureData;
    }
};

BLOCK.stoneBrick = {
    name: 'stone brick',
    id: 8,
    type: "cubeMesh",
    textureData: () =>
    {
        let textureData = [];
        textureData.push({texture: 'textures/stone_bricks.png', count: S.VERTEX_SIZE * 5});
        return textureData;
    }
};

BLOCK.glass = {
    name: 'glass',
    id: 9,
    type: "cubeMesh",
    textureData: () =>
    {
        let textureData = [];
        textureData.push({texture: 'textures/glass.png', count: S.VERTEX_SIZE * 5});
        return textureData;
    }
};

BLOCK.orchid = {
    name: 'flower',
    id: 10,
    type: "plantMesh",
    textureData: () =>
    {
        let textureData = [];
        textureData.push({texture: 'textures/blue_orchid.png', count: S.VERTEX_SIZE * 2});
        return textureData;
    }
};

BLOCK.flowerLeaf = {
    name: 'leaf block',
    id: 11,
    type: "cubeMesh",
    textureData: () =>
    {
        let textureData = [];
        textureData.push({texture: 'textures/flowering_azalea_leaves.png', count: S.VERTEX_SIZE * 5});
        return textureData;
    }
};

BLOCKFUNC.getBlock = (id) => {
    for ( let b in BLOCK )
        if ( BLOCK[b].id == id)
            return BLOCK[b];
    return null;
}

BLOCKFUNC.getBlockName = (id) => {
    for ( let b in BLOCK )
        if ( BLOCK[b].id == id) {
            return BLOCK[b].name;
        }
    return null;
}