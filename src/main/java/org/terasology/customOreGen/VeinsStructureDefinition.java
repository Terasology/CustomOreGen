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
import org.joml.Vector3i;
import org.joml.Vector3ic;
import org.joml.Vector4f;
import org.terasology.engine.utilities.random.FastRandom;
import org.terasology.engine.utilities.random.Random;
import org.terasology.joml.geom.AABBf;

import java.util.List;

/**
 * Code very heavily based on JRoush's implementation of CustomOreGen.
 * http://www.minecraftforum.net/topic/1107057-146v2-custom-ore-generation-updated-jan-5th/
 */
public class VeinsStructureDefinition extends AbstractMultiChunkStructureDefinition {

    private PDist motherLodeRadius;
    private PDist motherLodeYLevel;

    private PDist branchFrequency;
    private PDist branchInclination;
    private PDist branchLength;
    private PDist branchHeightLimit;

    private PDist segmentForkFrequency;
    private PDist segmentForkLengthMultiplier;
    private PDist segmentLength;
    private PDist segmentAngle;
    private PDist segmentRadius;

    private PDist blockDensity;
    private PDist blockRadiusMultiplier;

    public VeinsStructureDefinition(PDist frequency,
                                    PDist motherLodeRadius, PDist motherLodeYLevel,
                                    PDist branchFrequency, PDist branchInclination, PDist branchLength, PDist branchHeightLimit,
                                    PDist segmentForkFrequency, PDist segmentForkLengthMultiplier, PDist segmentLength, PDist segmentAngle, PDist segmentRadius,
                                    PDist blockDensity, PDist blockRadiusMultiplier) {
        super(frequency);
        this.motherLodeRadius = motherLodeRadius;
        this.motherLodeYLevel = motherLodeYLevel;
        this.branchFrequency = branchFrequency;
        this.branchInclination = branchInclination;
        this.branchLength = branchLength;
        this.branchHeightLimit = branchHeightLimit;
        this.segmentForkFrequency = segmentForkFrequency;
        this.segmentForkLengthMultiplier = segmentForkLengthMultiplier;
        this.segmentLength = segmentLength;
        this.segmentAngle = segmentAngle;
        this.segmentRadius = segmentRadius;
        this.blockDensity = blockDensity;
        this.blockRadiusMultiplier = blockRadiusMultiplier;
    }

    @Override
    protected int getGeneratorSalt() {
        return 2913452;
    }

    @Override
    protected float getMaxRange() {
        return motherLodeRadius.getMax() + branchLength.getMax();
    }

    @Override
    protected void generateStructuresForChunk(List<Structure> result, Random random, Vector3ic chunkSize, int xShift, int yShift, int zShift) {
        // motherlode X,Y,Z coordinates within chunk
        float mlX = random.nextFloat() * chunkSize.x() + xShift;
        float mlY = motherLodeYLevel.getValue(random) + yShift;
        float mlZ = random.nextFloat() * chunkSize.z() + zShift;

        // motherlode transformation matrix
        Matrix4f mlMat = new Matrix4f();
        mlMat.translate(mlX, mlY, mlZ); // center translation
        mlMat.rotateZ(random.nextFloat() * 6.28319F); // phi rotation
        mlMat.rotateY(random.nextFloat() * 6.28319F); // theta rotation
        mlMat.scale(motherLodeRadius.getValue(random), motherLodeRadius.getValue(random), motherLodeRadius.getValue(random)); // scale axes

        // create motherlode component
        result.add(new SolidSphereStructure(mlMat, chunkSize, random));

        // create random number of branches
        for (int br = branchFrequency.getIntValue(random); br > 0; br--) {
            // generate an independent random for this branch
            Random brRandom = new FastRandom(random.nextLong());
            // initialize branch transform
            Matrix4f segMat = new Matrix4f();
            segMat.translate(mlX, mlY, mlZ);  // motherlode translation
            segMat.rotateY(brRandom.nextFloat() * 6.28319F); // random rotation about vertical
            segMat.rotateX(-branchInclination.getValue(brRandom)); // angle from horizontal
            // calculate height limits
            float maxHeight = mlY + branchHeightLimit.getValue(brRandom);
            float minHeight = mlY - branchHeightLimit.getValue(brRandom);
            // create branch
            generateBranch(result, branchLength.getValue(brRandom), maxHeight, minHeight, segMat, null, brRandom);
        }
    }

    public void generateBranch(List<Structure> result, float length, float maxHeight, float minHeight, Matrix4f mat, BezierTubeStructure parent, Random random) {
        Vector4f pos = new Vector4f(0,0,0,0);
        float remainingLength = length;
        // create segments until max branch length is reached
        while (remainingLength > 0) {
            // determine segment length & radius
            float segLen = segmentLength.getValue(random);
            if (segLen > remainingLength) {
                segLen = remainingLength;
            }
            remainingLength -= segLen;
            segLen /= 2;
            float segRad = segmentRadius.getValue(random);

            // translate to center point
            mat.translate(0, 0, segLen);

            // create segment component
            Matrix4f segMat = new Matrix4f(mat).scale(segRad, segRad, segLen);
            BezierTubeStructure tube = new BezierTubeStructure(parent, segMat, random);
            result.add(tube);

            // translate to end point
            mat.translate(0, 0, segLen);
            // calculate end point
            pos.set(0,0,0,1);
            mat.transform(pos);

            // validate coordinates for next segment
            if (pos.y > maxHeight || pos.y < minHeight) {
                return;    // branch extends outside of vertical range
            } else if (remainingLength <= 0) {
                return; // remaining length is  <= 0
            }

            // create forks
            for (int fk = segmentForkFrequency.getIntValue(random); fk > 0; fk--) {
                // generate an independent random and transform for this fork
                Random fkRandom = new FastRandom(random.nextLong());
                Matrix4f fkMat = new Matrix4f(mat);
                // rotate relative to arbitrary axis in XY plane
                float axisTheta = fkRandom.nextFloat() * 6.28319F;
                fkMat.rotate(segmentAngle.getValue(fkRandom), (float) Math.cos(axisTheta), (float) Math.sin(axisTheta), 0);
                // create forked branch
                float fkLenMult = segmentForkLengthMultiplier.getValue(fkRandom);
                generateBranch(result, remainingLength * (fkLenMult > 1F ? 1F : fkLenMult), maxHeight, minHeight, fkMat, tube, fkRandom);
            }

            // rotate relative to arbitrary axis in XY plane
            float axisTheta = random.nextFloat() * 6.28319F;
            mat.rotate(segmentAngle.getValue(random), (float) Math.cos(axisTheta), (float) Math.sin(axisTheta), 0);
        }
    }

    private final class BezierTubeStructure implements Structure {
        // center and forward control points
        private Vector4f mid = new Vector4f(0f,0f,0f, 1.0f);
        private Vector4f end = new Vector4f(0f,0f,1f, 1.0f);
        // radius
        private final float rad;
        // neighbors
        private BezierTubeStructure prev;
        private Random random;
        private BezierTubeStructure next;
        // interpolation context & persistent transform object
        private final InterpolationContext context;
        private final Matrix4f mat;
        private Vector3i minPosition;
        private Vector3i maxPosition;

        private BezierTubeStructure(BezierTubeStructure parent, Matrix4f transform, Random random) {
            prev = parent;
            this.random = random;
            if (prev != null) {
                prev.next = this;
            }
            // calculate midpoint & endpoint
//            mid = new float[]{0, 0, 0};
            transform.transform(mid);
//            end = new float[]{0, 0, 1};
            transform.transform(end);
            // calculate radius (along x axis, we assume it is the same along the y)
            Vector4f xunit = new Vector4f(1, 0, 0, 0);
            transform.transform(xunit);
            rad = (float) Math.sqrt(xunit.x * xunit.x + xunit.y * xunit.y + xunit.z * xunit.z);
            // build transformed bounding box from the local BB for a cylinder
            float rMax = rad * blockRadiusMultiplier.getMax();
            if (rMax < 0) {
                rMax = 0;
            }
            AABBf bb = new AABBf(-rMax, -rMax, -1, rMax, rMax, 1);
            bb.transform(transform);



            minPosition = new Vector3i(Math.roundUsing(bb.minX, RoundingMode.FLOOR), Math.roundUsing(bb.minY,
                RoundingMode.FLOOR), Math.roundUsing(bb.minZ, RoundingMode.FLOOR));


            maxPosition = new Vector3i(Math.roundUsing(bb.maxX + 1, RoundingMode.FLOOR), Math.roundUsing(bb.maxY + 1,
                RoundingMode.FLOOR), Math.roundUsing(bb.maxZ + 1, RoundingMode.FLOOR));

            // construct a persistent context for interpolation loops
            context = new InterpolationContext();
            mat = new Matrix4f();
        }

        /**
         * Parametric interpolation of tube center line.
         * Quadratic bezier curves are used between neighboring segments.
         * Without a neighbor, the tube follows the straight segment line.
         *
         * @param pos set to position vector as {x,y,z}
         * @param t   interpolating parameter between [-1,1]
         */
        public void interpolatePosition(float[] pos, float t) {
            if (t > 0 && next != null) { // valid forward neighbor
                float nt = 1 - t;
                pos[0] = nt * nt * mid.x + 2 * t * nt * end.x + t * t * next.mid.x;
                pos[1] = nt * nt * mid.y + 2 * t * nt * end.y + t * t * next.mid.y;
                pos[2] = nt * nt * mid.z + 2 * t * nt * end.z + t * t * next.mid.z;
            } else if (t < 0 && prev != null) { // valid backward neighbor
                float nt = 1 + t;
                pos[0] = nt * nt * mid.x - 2 * t * nt * prev.end.x + t * t * prev.mid.x;
                pos[1] = nt * nt * mid.y - 2 * t * nt * prev.end.y + t * t * prev.mid.y;
                pos[2] = nt * nt * mid.z - 2 * t * nt * prev.end.z + t * t * prev.mid.z;
            } else { // no neighbor in specified direction - simple linear interpolation
                float nt = 1 - 2 * t;
                pos[0] = nt * mid.x + 2 * t * end.x;
                pos[1] = nt * mid.y + 2 * t * end.y;
                pos[2] = nt * mid.z + 2 * t * end.z;
            }
        }

        /**
         * Parametric interpolation of position derivative.
         * Without a neighbor, the tube follows the straight segment line.
         *
         * @param der set to derivative vector as {dx,dy,dz}
         * @param t   interpolating parameter between [-1,1]
         */
        public void interpolateDerivative(float[] der, float t) {
            if (t > 0 && next != null) { // valid forward neighbor
                der[0] = 2 * ((1 - t) * (end.x - mid.x) + t * (next.mid.x - end.x));
                der[1] = 2 * ((1 - t) * (end.y - mid.y) + t * (next.mid.y - end.y));
                der[2] = 2 * ((1 - t) * (end.z - mid.z) + t * (next.mid.z - end.z));
            } else if (t < 0 && prev != null) { // valid backward neighbor
                der[0] = 2 * ((1 + t) * (mid.x - prev.end.x) - t * (prev.end.x - prev.mid.x));
                der[1] = 2 * ((1 + t) * (mid.y - prev.end.y) - t * (prev.end.y - prev.mid.y));
                der[2] = 2 * ((1 + t) * (mid.z - prev.end.z) - t * (prev.end.z - prev.mid.z));
            } else { // no neighbor in specified direction
                der[0] = 2 * (end.x - mid.x);
                der[1] = 2 * (end.y - mid.y);
                der[2] = 2 * (end.z - mid.z);
            }
        }

        /**
         * Parametric interpolation of segment radius.
         * Segment radius varies smoothly between neighboring segments.
         * Without a neighbor, the forward half terminates in a half-ellipsoid
         * and the backward half terminates in a cylinder of constant radius
         *
         * @param t interpolating parameter between [-1,1]
         */
        public float interpolateRadius(float t) {
            if (t > 0 && next != null) {
                return (1 - t) * rad + t * next.rad; // valid forward neighbor
            } else if (t < 0 && prev != null) {
                return (1 + t) * rad - t * prev.rad; // valid backward neighbor
            } else if (t <= 0 && t > -1) {
                return rad; // no backward neighbor - constant radius
            } else if (t > 0 && t < 1) {
                return rad * (float) Math.sqrt(1 - 4 * t * t); // no forward neighbor - approach zero as parabola
            } else {
                return 0;
            }
        }

        @Override
        public void generateStructure(StructureCallback callback) {
            // get min & max radii in local coordinates
            float maxR = blockRadiusMultiplier.getMax();
            if (maxR < 0) {
                maxR = 0;
            }
            float maxR2 = maxR * maxR;
            float minR = blockRadiusMultiplier.getMin();
            if (minR < 0) {
                minR = 0;
            }
            float minR2 = minR * minR;

            // interpolate over segment
            Vector4f pos = new Vector4f();
            int innerStep = 1;
            context.init(0, true);
            do {
                // determine step size and count
                innerStep = (int) context.radius / 4 + 1;
                if (context.radius <= 0) {
                    continue; // zero radius
                }
                float step = 0.7F * innerStep / context.radius;
                int stepCount = (int) (maxR / step) + 1;
                boolean oneBlockThreshold = (context.radius * maxR < 0.25F); // radius is too small even for a single block
                // build transformation
                mat.identity();
                mat.translate(context.pos[0], context.pos[1], context.pos[2]);
                mat.rotateTowards(0,1,0, context.der[0], context.der[1], context.der[2]);
                mat.scale(context.radius, context.radius, innerStep);

                Vector3i blockPos = new Vector3i();
                Vector3i blockCenter = new Vector3i();
                // iterate through blocks in the local XY plane
                for (int x = -stepCount; x < stepCount; x++) {
                    for (int y = -stepCount; y < stepCount; y++) {
                        pos.x = x * step;
                        pos.y = y * step;
                        pos.z = 0;
                        // check radius
                        float r2 = pos.x * pos.x + pos.y * pos.y;
                        if (r2 > maxR2) {
                            continue; // block is outside maximum possible radius
                        }
                        if (r2 > minR2) { // block is near tube surface
                            float rMax = blockRadiusMultiplier.getValue(random);
                            if (r2 > rMax * rMax) {
                                continue; // block is outside maximum radius
                            }
                        }
                        if (oneBlockThreshold && context.radius * maxR * 4 < random.nextFloat()) {
                            continue; // blocks must pass random check for very thin tubes
                        }
                        // transform into world coordinates
                        mat.transform(pos);
                        int baseX = (int) Math.floor(pos.x) - innerStep / 2;
                        int baseY = (int) Math.floor(pos.y) - innerStep / 2;
                        int baseZ = (int) Math.floor(pos.z) - innerStep / 2;
                        // iterate over inner group
                        for (int blockX = baseX; blockX < innerStep + baseX; blockX++) {
                            for (int blockY = baseY; blockY < innerStep + baseY; blockY++) {
                                for (int blockZ = baseZ; blockZ < innerStep + baseZ; blockZ++) {
                                    if (callback.canReplace(blockX, blockY, blockZ) && blockDensity.getIntValue(random) >= 1) {
                                        callback.replaceBlock(blockPos.set(blockX, blockY, blockZ), StructureNodeType.BRANCH, blockCenter);
                                    }
                                }
                            }
                        }
                    }
                }
                // next interpolation step
            } while (context.advance(0.7F * innerStep));
        }

        /**
         * Context information for discrete interpolation over the segement
         */
        private final class InterpolationContext {
            public float[] pos;        // position
            public float[] der;        // normalized derivative vector
            public float derLen;    // norm of derivative
            public float radius;    // radius of tube
            public float err;        // estimated max distance^2 to corresponding points from previous step
            public float t;            // interpolation parameter
            public float dt;        // parameter step size
            public boolean calcDer;    // whether or not to calculate the derivative vector and length at each step

            /**
             * Construct a blank context object
             */
            public InterpolationContext() {
                pos = new float[3];
                der = new float[3];
                t = 10; // interpolation is not in progress
                dt = 1 / 20F; // guess at step size
            }

            /**
             * Called before an interpolation loop to initialize all state info
             *
             * @param stepSize           Initial step size to use.  Pass 0 to keep size from last interpolation.
             * @param calculateDirection Whether or not to compute normalized derivative at each step.
             */
            public void init(float stepSize, boolean calculateDirection) {
                // interpolate all the way to the center of previous segment unless it reciprocates
                t = (prev == null || prev.next == BezierTubeStructure.this) ? -0.5F : -1.0F;
                if (stepSize > 0) {
                    dt = stepSize;
                }
                // calculate initial position
                interpolatePosition(pos, t);
                // calculate initial radius
                radius = interpolateRadius(t);
                // calculate and normalize initial derivative
                calcDer = calculateDirection;
                if (calcDer) {
                    interpolateDerivative(der, t);
                    derLen = (float) Math.sqrt(der[0] * der[0] + der[1] * der[1] + der[2] * der[2]);
                    der[0] /= derLen;
                    der[1] /= derLen;
                    der[2] /= derLen;
                } else {
                    derLen = 0;
                    der[0] = 0;
                    der[1] = 0;
                    der[2] = 0;
                }
                // set initial error to zero
                err = 0;
            }

            /**
             * Called to advance a step in the interpolation
             *
             * @return Returns false when the interpolation has reached the end of the segment
             */
            public boolean advance(float tolerance) {
                // store current state
                float pX = pos[0];
                float pY = pos[1];
                float pZ = pos[2];
                float dX = der[0];
                float dY = der[1];
                float dZ = der[2];
                float r = radius;
                // attempt to advance
                while (true) {
                    // advance parameter
                    float nt = t + dt;
                    // calculate new position and displacement^2
                    interpolatePosition(pos, nt);
                    float deltaX = pX - pos[0];
                    float deltaY = pY - pos[1];
                    float deltaZ = pZ - pos[2];
                    float d2 = deltaX * deltaX + deltaY * deltaY + deltaZ * deltaZ;
                    err = d2; // displacement ^2 due to translation along curve
                    // calculate new radius
                    radius = interpolateRadius(nt);
                    float avg2R = r + radius;
                    // calculate new derivative and sin^2 of angle between old & new derivatives
                    if (calcDer) {
                        interpolateDerivative(der, nt);
                        derLen = (float) Math.sqrt(der[0] * der[0] + der[1] * der[1] + der[2] * der[2]);
                        der[0] /= derLen;
                        der[1] /= derLen;
                        der[2] /= derLen;
                        deltaX = -dZ * der[1] + dY * der[2];
                        deltaY = dZ * der[0] - dX * der[2];
                        deltaZ = -dY * der[0] + dX * der[1];
                        float sin2 = deltaX * deltaX + deltaY * deltaY + deltaZ * deltaZ;
                        err += avg2R * avg2R * sin2; // ~ max displacement^2 due to relative rotation
                    }
                    // check error
                    float maxErr = tolerance * tolerance;
                    if (err > maxErr) {
                        dt *= 0.6;  // reduce step size -> reduce error
                    } else if (err < maxErr / 5) {
                        dt *= 1.8;    // increase step size -> fewer steps
                    } else {
                        break; // error was acceptable
                    }
                    // prevent infinite loops
                    if (dt < java.lang.Math.ulp(t) * 2) {
                        throw new RuntimeException("CustomOreGen: Detected a possible infinite loop during bezier interpolation.  Please report this error.");
                    }
                }
                // advancement succeeded - advance actual parameter
                t += dt;
                return (t < 0.5F);
            }
        }
    }


    private final class SolidSphereStructure implements Structure {
        protected final Matrix4f mat;
        protected final Matrix4f invMat;
        private Vector3i minPosition;
        private Vector3i maxPosition;
        private Vector3ic chunkSize;
        private Random random;

        public SolidSphereStructure(Matrix4f transform, Vector3ic chunkSize, Random random) {
            this.chunkSize = chunkSize;
            this.random = random;
            // build transformed bounding box from the local BB for a unit sphere
            float rMax = blockRadiusMultiplier.getMax();
            if (rMax < 0) {
                rMax = 0;
            }

            AABBf aabb = new AABBf(-rMax, -rMax, -rMax, rMax, rMax, rMax);
            aabb.transform(transform);
            minPosition = new Vector3i(Math.roundUsing(aabb.minX, RoundingMode.FLOOR), Math.roundUsing(aabb.minY, RoundingMode.FLOOR), Math.roundUsing(aabb.minZ, RoundingMode.FLOOR));
            maxPosition = new Vector3i(Math.roundUsing(aabb.maxX + 1, RoundingMode.FLOOR), Math.roundUsing(aabb.maxY + 1, RoundingMode.FLOOR), Math.roundUsing(aabb.maxZ + 1, RoundingMode.FLOOR));

            // store transforms
            mat = new Matrix4f(transform);
            if (transform.determinant() != 0) {
                invMat = transform.invert(new Matrix4f());
            } else {
                invMat = null; // at least one axis of sphere has zero length
            }
        }

        @Override
        public void generateStructure(StructureCallback callback) {
            if (invMat == null) {
                return; // sphere has zero volume and therefore cannot contain blocks
            }
            // get min & max radii in local coordinates
            float maxR2 = blockRadiusMultiplier.getMax();
            if (maxR2 < 0) {
                maxR2 = 0;
            }
            maxR2 *= maxR2;
            float minR2 = blockRadiusMultiplier.getMin();
            if (minR2 < 0) {
                minR2 = 0;
            }
            minR2 *= minR2;
            // iterate through blocks
            Vector4f pos = new Vector4f();
            for (int x = Math.max(0, minPosition.x); x <= Math.min(chunkSize.x() - 1, maxPosition.x); x++) {
                for (int y = Math.max(0, minPosition.y); y <= Math.min(chunkSize.y() - 1, maxPosition.y); y++) {
                    for (int z = Math.max(0, minPosition.z); z <= Math.min(chunkSize.z() - 1, maxPosition.z); z++) {
                        if (!callback.canReplace(x, y, z)) {
                            continue;
                        }
                        // transform into local coordinates
                        pos.x = x + 0.5F;
                        pos.y = y + 0.5F;
                        pos.z = z + 0.5F;
                        invMat.transform(pos);
                        // check radius
                        float r2 = pos.x * pos.x + pos.y * pos.y + pos.z * pos.z;
                        if (r2 > maxR2) {
                            continue; // block is outside maximum possible radius
                        }
                        if (r2 > minR2) { // block is near ellipsoid surface
                            float rMax = blockRadiusMultiplier.getValue(random);
                            if (r2 > rMax * rMax) {
                                continue; // block is outside maximum radius
                            }
                        }
                        // place block
                        if (blockDensity.getIntValue(random) < 1) {
                            continue; // density check failed
                        }

                        Vector3i blockPosition = new Vector3i(x, y, z);
                        callback.replaceBlock(blockPosition, StructureNodeType.CLUSTER, getRelativePosition(blockPosition, minPosition, maxPosition));
                    }
                }
            }
        }
    }
}
