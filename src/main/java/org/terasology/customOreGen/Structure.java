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

import org.joml.Vector3ic;

public interface Structure {
    void generateStructure(StructureCallback callback);

    public interface StructureCallback {
        /**
         * Relative position
         */
        void replaceBlock(Vector3ic position, StructureNodeType structureNodeType, Vector3ic distanceToCenter);

        /**
         * Relative position
         *
         * @param x
         * @param y
         * @param z
         * @return
         */
        boolean canReplace(int x, int y, int z);
    }
}
