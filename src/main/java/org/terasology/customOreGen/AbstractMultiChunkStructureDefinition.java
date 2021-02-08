// Copyright 2021 The Terasology Foundation
// SPDX-License-Identifier: Apache-2.0
package org.terasology.customOreGen;

import org.joml.Vector3f;
import org.joml.Vector3i;
import org.joml.Vector3ic;
import org.terasology.utilities.random.Random;
import org.terasology.world.block.BlockRegion;
import org.terasology.world.chunks.Chunks;

import java.util.Collection;
import java.util.LinkedList;
import java.util.List;

public abstract class AbstractMultiChunkStructureDefinition implements StructureDefinition {
    private PDist frequency;

    protected AbstractMultiChunkStructureDefinition(PDist frequency) {
        this.frequency = frequency;
    }

    @Override
    public final Collection<Structure> generateStructures(long seed, BlockRegion worldRegion) {
        List<Structure> result = new LinkedList<>();
        Vector3i chunkPosition = Chunks.toChunkPos(worldRegion.center(new Vector3f()), new Vector3i());
        float maxRange = getMaxRange();
        Vector3i chunkSize = worldRegion.getSize(new Vector3i());

        int horzontalChunkRange = (int) Math.ceil(Math.max(maxRange / chunkSize.x, maxRange / chunkSize.z));
        int verticalChunkRange = (int) Math.ceil(maxRange / chunkSize.y);
        for (int chunkX = -horzontalChunkRange; chunkX <= horzontalChunkRange; chunkX++) {
            for (int chunkZ = -horzontalChunkRange; chunkZ <= horzontalChunkRange; chunkZ++) {
                for (int chunkY = -verticalChunkRange; chunkY <= verticalChunkRange; chunkY++) {

                    generateStructuresForChunkWithFrequency(result, seed,
                        new Vector3i(chunkPosition.x + chunkX, chunkPosition.y + chunkY, chunkPosition.z + chunkZ),
                        chunkSize, chunkX * chunkSize.x, chunkY * chunkSize.y, chunkZ * chunkSize.z);
                }
            }
        }

        return result;
    }

    protected final void generateStructuresForChunkWithFrequency(List<Structure> result, long seed, Vector3ic chunkPosition,
                                                                 Vector3ic chunkSize, int xShift, int yShift, int zShift) {
        Random random = ChunkRandom.getChunkRandom(seed, chunkPosition, getGeneratorSalt());

        float structuresInChunk = frequency.getValue(random);
        int structuresToGenerateInChunk = (int) structuresInChunk;

        // Check if we "hit" any leftover
        if (random.nextFloat() < structuresInChunk - structuresToGenerateInChunk) {
            structuresToGenerateInChunk++;
        }

        for (int i = 0; i < structuresToGenerateInChunk; i++) {
            generateStructuresForChunk(result, random, chunkSize, xShift, yShift, zShift);
        }
    }

    protected abstract float getMaxRange();

    protected abstract int getGeneratorSalt();

    protected abstract void generateStructuresForChunk(List<Structure> result, Random random, Vector3ic chunkSize, int xShift, int yShift, int zShift);

    protected Vector3i getRelativePosition(Vector3ic blockPosition, Vector3ic originMinPosition, Vector3ic originMaxPosition) {
        Vector3i relativePosition = new Vector3i(originMaxPosition);
        relativePosition.add(originMinPosition);
        relativePosition.div(2);
        relativePosition.sub(blockPosition);
        return relativePosition;
    }
}
