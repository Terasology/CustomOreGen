/*
 * Copyright 2014 MovingBlocks
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package org.terasology.customOreGen;

import org.terasology.math.geom.Vector3i;
import org.terasology.utilities.random.FastRandom;
import org.terasology.utilities.random.Random;

public abstract class ChunkRandom {
    public static Random getChunkRandom(long seed, Vector3i location, int salt) {
        return new FastRandom(seed + salt * (97 * location.x + location.y + location.z));
    }
}
