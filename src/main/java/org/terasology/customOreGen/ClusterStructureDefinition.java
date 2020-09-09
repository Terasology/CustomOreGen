// Copyright 2020 The Terasology Foundation
// SPDX-License-Identifier: Apache-2.0
package org.terasology.customOreGen;

import org.terasology.engine.utilities.random.Random;
import org.terasology.math.geom.Vector3i;

import java.util.List;

/**
 * Code very heavily based on JRoush's implementation of CustomOreGen. http://www.minecraftforum
 * .net/topic/1107057-146v2-custom-ore-generation-updated-jan-5th/
 */
public class ClusterStructureDefinition extends AbstractMultiChunkStructureDefinition {
    private final PDist pocketYLevel;
    private final PDist clusterRichness;

    public ClusterStructureDefinition(PDist frequency, PDist clusterRichness, PDist pocketYLevel) {
        super(frequency);
        this.clusterRichness = clusterRichness;
        this.pocketYLevel = pocketYLevel;
    }

    @Override
    protected int getGeneratorSalt() {
        return 979237;
    }

    @Override
    protected float getMaxRange() {
        return clusterRichness.getMax();
    }

    @Override
    protected void generateStructuresForChunk(List<Structure> result, Random random, Vector3i chunkSize, int xShift,
                                              int yShift, int zShift) {
        // cluster X,Y,Z coordinates within chunk
        float clX = random.nextFloat() * chunkSize.x + xShift;
        float clY = pocketYLevel.getValue(random) + yShift;
        float clZ = random.nextFloat() * chunkSize.z + zShift;

        result.add(new ClusterStructure(clX, clY, clZ, random, chunkSize));
    }

    private class ClusterStructure implements Structure {
        /* Cluster geometry parameters */
        protected final int size;        // scale
        protected final float[] ptA;    // segment start
        protected final float[] ptB;    // segment end
        protected final float[] rad;    // radius at each step along segment
        private final Vector3i chunkSize;

        public ClusterStructure(float x, float y, float z, Random random, Vector3i chunkSize) {
            this.chunkSize = chunkSize;
            // choose segment length and horizontal angle from +Z axis
            size = clusterRichness.getIntValue(random);
            double horizAngle = random.nextFloat() * Math.PI;
            ptA = new float[3];
            ptB = new float[3];

            // rotate segment in XZ plane
            float segmentXOffset = (float) Math.sin(horizAngle) * size / 8F;
            float segmentZOffset = (float) Math.cos(horizAngle) * size / 8F;
            ptA[0] = x + segmentXOffset;
            ptB[0] = x - segmentXOffset;
            ptA[2] = z + segmentZOffset;
            ptB[2] = z - segmentZOffset;

            // drop each endpoint randomly by up to 2 blocks
            // not sure why only negative changes are allowed, possibly
            // in deference to the upper height limits imposed on ore spawning.
            ptA[1] = y + random.nextInt(3) - 2;
            ptB[1] = y + random.nextInt(3) - 2;

            // compute radii along segment and expand BB
            // Radius has a random factor from [0,size/32] and a position factor that smoothly
            // increases the radius towards the center of the segment, up of a max of x2 at the
            // middle.
            rad = new float[size + 1];
            for (int s = 0; s < rad.length; s++) {
                // radius for this step
                float ns = s / (float) (rad.length - 1);
                float baseRadius = (float) random.nextDouble() * size / 32F;
                rad[s] = ((float) Math.sin(ns * Math.PI) + 1) * baseRadius + 0.5F;
            }
        }

        @Override
        public void generateStructure(StructureCallback callback) {
            // iterate along segment in unit steps
            for (int s = 0; s < rad.length; s++) {
                float ns = s / (float) (rad.length - 1);

                // center of sphere for this step
                float xCenter = ptA[0] + (ptB[0] - ptA[0]) * ns;
                float yCenter = ptA[1] + (ptB[1] - ptA[1]) * ns;
                float zCenter = ptA[2] + (ptB[2] - ptA[2]) * ns;

                // iterate over each block in the bounding box of the sphere
                int xMin = (int) Math.max(0, Math.floor(xCenter - rad[s]));
                int xMax = (int) Math.min(chunkSize.x - 1, Math.ceil(xCenter + rad[s]));
                int yMin = (int) Math.max(0, Math.floor(yCenter - rad[s]));
                int yMax = (int) Math.min(chunkSize.y - 1, Math.ceil(yCenter + rad[s]));
                int zMin = (int) Math.max(0, Math.floor(zCenter - rad[s]));
                int zMax = (int) Math.min(chunkSize.z - 1, Math.ceil(zCenter + rad[s]));
                for (int tgtX = xMin; tgtX <= xMax; tgtX++) {
                    double normXDist = (tgtX + 0.5D - xCenter) / rad[s];
                    if (normXDist * normXDist >= 1.0D) {
                        continue;
                    }
                    for (int tgtY = yMin; tgtY <= yMax; tgtY++) {
                        double normYDist = (tgtY + 0.5D - yCenter) / rad[s];
                        if (normXDist * normXDist + normYDist * normYDist >= 1.0D) {
                            continue;
                        }
                        for (int tgtZ = zMin; tgtZ <= zMax; tgtZ++) {
                            double normZDist = (tgtZ + 0.5D - zCenter) / rad[s];
                            if (normXDist * normXDist + normYDist * normYDist + normZDist * normZDist >= 1.0D) {
                                continue;
                            }

                            callback.replaceBlock(new Vector3i(tgtX, tgtY, tgtZ), StructureNodeType.CLUSTER,
                                    Vector3i.zero());
                        }
                    }
                }
            }
        }
    }
}


