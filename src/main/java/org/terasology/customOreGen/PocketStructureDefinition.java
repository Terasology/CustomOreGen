/*
 * Copyright 2014 MovingBlocks
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *  http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package org.terasology.customOreGen;

import org.joml.Math;
import org.joml.Matrix4f;
import org.joml.RoundingMode;
import org.joml.Vector3f;
import org.joml.Vector3i;
import org.joml.Vector3ic;
import org.joml.Vector4f;
import org.terasology.utilities.procedural.Noise3D;
import org.terasology.utilities.procedural.SimplexNoise;
import org.terasology.utilities.random.Random;
import org.terasology.world.block.Block;

import java.util.List;

/**
 * Code very heavily based on JRoush's implementation of CustomOreGen.
 * http://www.minecraftforum.net/topic/1107057-146v2-custom-ore-generation-updated-jan-5th/
 */
public class PocketStructureDefinition extends AbstractMultiChunkStructureDefinition {
    private PDist pocketRadius;
    private PDist pocketThickness;
    private PDist pocketYLevel;
    private PDist pocketAngle;
    private PDist blockRadiusMult;
    private PDist blockDensity;
    private PDist noiseLevel;
    private PDist volumeNoiseCutOff;

    public PocketStructureDefinition(PDist frequency, PDist pocketRadius, PDist pocketThickness, PDist pocketYLevel, PDist pocketAngle,
                                     PDist blockRadiusMult, PDist blockDensity, PDist noiseLevel, PDist volumeNoiseCutOff) {
        super(frequency);
        this.pocketRadius = pocketRadius;
        this.pocketThickness = pocketThickness;
        this.pocketYLevel = pocketYLevel;
        this.pocketAngle = pocketAngle;
        this.blockRadiusMult = blockRadiusMult;
        this.blockDensity = blockDensity;
        this.noiseLevel = noiseLevel;
        this.volumeNoiseCutOff = volumeNoiseCutOff;
    }

    protected float getMaxRange() {
        return Math.max(pocketRadius.getMax(), pocketThickness.getMax());
    }

    @Override
    protected int getGeneratorSalt() {
        return 3423876;
    }

    @Override
    protected void generateStructuresForChunk(List<Structure> result, Random random, Vector3ic chunkSize, int xShift, int yShift, int zShift) {
        // cloud X,Y,Z coordinates within chunk
        float clX = random.nextFloat() * chunkSize.x() + xShift;
        float clY = pocketYLevel.getValue(random) + yShift;
        float clZ = random.nextFloat() * chunkSize.z() + zShift;

        // cloud transformation matrix
        Matrix4f clMat = new Matrix4f();
        clMat.translation(clX, clY, clZ);
        clMat.rotateTowards(0, 1, 0, 0, 0, 1);
        clMat.rotateZ(random.nextFloat() * 6.28319F);
        clMat.rotateY(pocketAngle.getValue(random));
        clMat.scale(pocketRadius.getValue(random), pocketRadius.getValue(random), pocketThickness.getValue(random));

        // create cloud component
        result.add(new DiffusePocketStructure(clMat, random, chunkSize));
    }

    public interface PocketBlockProvider {
        Block getBlock(float distanceFromCenter);
    }

    /**
     * A diffuse cloud of ore.
     */
    private class DiffusePocketStructure implements Structure {
        // transformation matrices
        protected final Matrix4f mat;
        protected final Matrix4f invMat;
        // noise generator
        protected final Noise3D noiseGen;
        protected final float sizeNoiseMagnitude;
        protected final int noiseLevels;
        private Vector3ic minPosition;
        private Vector3ic maxPosition;
        private Random random;
        private Vector3ic chunkSize;

        public DiffusePocketStructure(Matrix4f transform, Random random, Vector3ic chunkSize) {
            this.random = random;
            this.chunkSize = chunkSize;
            // create noise generator
            noiseGen = new SimplexNoise(random.nextInt());
            sizeNoiseMagnitude = Math.abs(noiseLevel.getValue(random));

            // build transformed bounding box from the local BB for a unit sphere
            float rMax = (1 + sizeNoiseMagnitude * 2) * blockRadiusMult.getMax();
            if (rMax < 0) {
                rMax = 0;
            }
            Vector3f min = new Vector3f();
            Vector3f max = new Vector3f();
            transform.scale(rMax).frustumAabb(min, max);

            minPosition = new Vector3i(Math.roundUsing(min.x, RoundingMode.FLOOR), Math.roundUsing(min.y, RoundingMode.FLOOR), Math.roundUsing(min.z, RoundingMode.FLOOR));
            maxPosition = new Vector3i(Math.roundUsing(max.x + 1, RoundingMode.FLOOR), Math.roundUsing(max.y + 1, RoundingMode.FLOOR), Math.roundUsing(max.z + 1, RoundingMode.FLOOR));

            // calculate noise levels from size of BB
            float maxSize = Math.max(max.x - min.x, Math.max(max.y - min.y, max.z - min.z)) * 0.2F;
            noiseLevels = (maxSize <= 1) ? 0 : (int) (java.lang.Math.log(maxSize) / java.lang.Math.log(2) + 0.5F);

            // store transforms
            mat = new Matrix4f(transform);
            if (transform.determinant() != 0) {
                invMat = transform.invert(new Matrix4f());
            } else {
                invMat = null; // at least one axis of sphere has zero length
            }
        }

        /**
         * Get total 1/f noise value at the specified position
         */
        public float getNoise(float x, float y, float z) {
            double noise = 0;
            for (int i = 0; i < noiseLevels; i++) {
                float im = (1 << i);
                noise += (1 / im) * noiseGen.noise(x * im, y * im, z * im); // add 1/f noise
            }
            return (float) noise;
        }

        @Override
        public void generateStructure(StructureCallback callback) {
            if (invMat == null) {
                return; // sphere has zero volume and therefore cannot contain blocks
            }

            // get min & max radii in local coordinates
            float maxR = Math.max(blockRadiusMult.getMax(), 0); // maximum radius after noise scaling
            float minR = Math.max(blockRadiusMult.getMin(), 0); // minimum radius after noise scaling
            float maxNoisyR2 = maxR * (1 + sizeNoiseMagnitude * 2); // maximum radius before noise scaling
            float minNoisyR2 = minR * (1 - sizeNoiseMagnitude * 2); // minimum radius before noise scaling
            maxNoisyR2 *= maxNoisyR2;
            minNoisyR2 *= minNoisyR2;
            // iterate through blocks
            Vector4f pos = new Vector4f();
            for (int x = Math.max(0, minPosition.x()); x <= Math.min(chunkSize.x() - 1, maxPosition.x()); x++) {
                for (int y = Math.max(0, minPosition.y()); y <= Math.min(chunkSize.y() - 1, maxPosition.y()); y++) {
                    for (int z = Math.max(0, minPosition.z()); z <= Math.min(chunkSize.z() - 1, maxPosition.z()); z++) {
                        if (!callback.canReplace(x, y, z)) {
                            continue;
                        }

                        // transform into local coordinates
                        pos.x = x + 0.5F;
                        pos.y = y + 0.5F;
                        pos.z = z + 0.5F;
                        invMat.transform(pos);
                        // check radius
                        float r2 = pos.lengthSquared();
                        if (r2 > maxNoisyR2) {
                            continue; // block is outside maximum possible radius
                        }
                        if (r2 > minNoisyR2) { // block is within max noise tolerance
                            // compute radius noise multiplier
                            // the point is projected radially onto the surface of the unit sphere and
                            // the 3D noise field is sampled at that location to determine the radius
                            // multiplier for that solid angle.
                            float r = (float) Math.sqrt(r2);
                            float mult = 1;
                            if (r > 0) {
                                mult += sizeNoiseMagnitude * getNoise(pos.x / r, pos.y / r, pos.z / r);
                            } else {
                                mult += sizeNoiseMagnitude * getNoise(0, 0, 0);
                            }
                            if (mult <= 0) {
                                continue; // noise-multiplied radius at this solid angle is zero
                            }
                            r /= mult;
                            // check noise-scaled radius
                            if (r > maxR) {
                                continue; // unit radius is outside max cutoff
                            }
                            if (r > minR && r > blockRadiusMult.getValue(random)) {
                                continue; // block is outside cutoff
                            }
                        }
                        // apply internal density noise
                        if (volumeNoiseCutOff.getMin() > 1) {
                            continue; // noise cutoff is too high
                        } else if (volumeNoiseCutOff.getMax() > 0) {
                            if ((getNoise(pos.x, pos.y, pos.z) + 1) / 2 < volumeNoiseCutOff.getValue(random)) {
                                continue; // noise level below cutoff
                            }
                        }
                        if (blockDensity.getIntValue(random) < 1) {
                            continue; // density check failed
                        }

                        // place block
                        Vector3i blockPosition = new Vector3i(x, y, z);
                        callback.replaceBlock(blockPosition, StructureNodeType.POCKET, getRelativePosition(blockPosition, minPosition, maxPosition));
                    }
                }
            }
        }
    }
}
