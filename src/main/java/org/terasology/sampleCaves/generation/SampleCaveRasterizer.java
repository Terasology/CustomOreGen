// Copyright 2020 The Terasology Foundation
// SPDX-License-Identifier: Apache-2.0
package org.terasology.sampleCaves.generation;

import org.terasology.engine.registry.CoreRegistry;
import org.terasology.engine.world.block.Block;
import org.terasology.engine.world.block.BlockManager;
import org.terasology.engine.world.chunks.ChunkConstants;
import org.terasology.engine.world.chunks.CoreChunk;
import org.terasology.engine.world.generation.Region;
import org.terasology.engine.world.generation.WorldRasterizer;
import org.terasology.math.geom.Vector3i;

public class SampleCaveRasterizer implements WorldRasterizer {
    String blockUri;

    public SampleCaveRasterizer() {
    }

    public SampleCaveRasterizer(String blockUri) {
        this.blockUri = blockUri;
    }

    @Override
    public void initialize() {
    }

    @Override
    public void generateChunk(CoreChunk chunk, Region chunkRegion) {
        SampleCaveFacet sampleCaveFacet = chunkRegion.getFacet(SampleCaveFacet.class);

        BlockManager blockManager = CoreRegistry.get(BlockManager.class);
        Block caveBlock = blockManager.getBlock(BlockManager.AIR_ID);
        if (blockUri != null) {
            caveBlock = blockManager.getBlock(blockUri);
        }

        for (Vector3i position : ChunkConstants.CHUNK_REGION) {
            if (sampleCaveFacet.get(position)) {
                chunk.setBlock(position, caveBlock);
            }
        }
    }
}
