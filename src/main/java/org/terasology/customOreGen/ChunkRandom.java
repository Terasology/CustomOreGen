// Copyright 2020 The Terasology Foundation
// SPDX-License-Identifier: Apache-2.0
package org.terasology.customOreGen;

import org.terasology.engine.utilities.random.FastRandom;
import org.terasology.engine.utilities.random.Random;
import org.terasology.math.geom.Vector3i;

public abstract class ChunkRandom {
    public static Random getChunkRandom(long seed, Vector3i location, int salt) {
        return new FastRandom(seed + salt * (97 * location.x + location.y + location.z));
    }
}
