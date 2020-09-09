// Copyright 2020 The Terasology Foundation
// SPDX-License-Identifier: Apache-2.0
package org.terasology.customOreGen;

import org.terasology.engine.math.Region3i;

import java.util.Collection;

public interface StructureDefinition {
    Collection<Structure> generateStructures(long seed, Region3i worldRegion);
}
