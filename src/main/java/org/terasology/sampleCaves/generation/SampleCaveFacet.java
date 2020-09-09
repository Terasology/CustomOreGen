// Copyright 2020 The Terasology Foundation
// SPDX-License-Identifier: Apache-2.0
package org.terasology.sampleCaves.generation;

import org.terasology.engine.math.Region3i;
import org.terasology.engine.world.generation.Border3D;
import org.terasology.engine.world.generation.facets.base.BaseBooleanFieldFacet3D;

public class SampleCaveFacet extends BaseBooleanFieldFacet3D {
    public SampleCaveFacet(Region3i targetRegion, Border3D border) {
        super(targetRegion, border);
    }
}
