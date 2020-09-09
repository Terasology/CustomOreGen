// Copyright 2020 The Terasology Foundation
// SPDX-License-Identifier: Apache-2.0
package org.terasology.sampleCaves.generation;

import org.terasology.customOreGen.Structure;
import org.terasology.customOreGen.StructureNodeType;
import org.terasology.engine.world.generation.facets.base.BaseBooleanFieldFacet3D;
import org.terasology.math.geom.Vector3i;

public class SampleCaveStructureCallback implements Structure.StructureCallback {
    private final BaseBooleanFieldFacet3D facet;

    public SampleCaveStructureCallback(BaseBooleanFieldFacet3D facet) {
        this.facet = facet;
    }

    @Override
    public boolean canReplace(int x, int y, int z) {
        return facet.getRelativeRegion().encompasses(x, y, z);
    }

    @Override
    public void replaceBlock(Vector3i position, StructureNodeType structureNodeType, Vector3i distanceToCenter) {
        if (canReplace(position.x, position.y, position.z)) {
            facet.set(position, true);
        }
    }
}
