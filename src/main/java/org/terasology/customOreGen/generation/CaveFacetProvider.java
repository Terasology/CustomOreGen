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
package org.terasology.customOreGen.generation;

import org.terasology.customOreGen.PDist;
import org.terasology.customOreGen.Structure;
import org.terasology.customOreGen.StructureDefinition;
import org.terasology.customOreGen.VeinsStructureDefinition;
import org.terasology.world.generation.FacetProvider;
import org.terasology.world.generation.GeneratingRegion;
import org.terasology.world.generation.Produces;

import java.util.Collection;

@Produces(CaveFacet.class)
public class CaveFacetProvider implements FacetProvider {
    long seed;

    @Override
    public void setSeed(long seed) {

    }

    @Override
    public void process(GeneratingRegion region) {
        CaveFacet facet = new CaveFacet(region.getRegion(), region.getBorderForFacet(CaveFacet.class));

        /* AnotherWorld
        PDist caveFrequency = new PDist(0.1f, 0f);
        PDist mainCaveRadius = new PDist(5f, 1f);
        PDist mainCaveYLevel = new PDist(32f, 32f);
        PDist tunnelLength = new PDist(45f, 10f);
        PDist tunnelRadius = new PDist(2f, 0.5f);
        */

        // note, try to contain the maxRange to an optimal range of neighboring chunks
        PDist caveFrequency = new PDist(0.1f, 0f);
        PDist mainCaveRadius = new PDist(8f, 3f);
        PDist mainCaveYLevel = new PDist(32f, 32f);
        PDist tunnelLength = new PDist(40f, 10f);
        PDist tunnelRadius = new PDist(4f, 1f);

        StructureDefinition structureDefinition = new VeinsStructureDefinition(caveFrequency,
                mainCaveRadius, mainCaveYLevel, new PDist(4f, 1f), new PDist(0f, 0.1f), tunnelLength,
                new PDist(100f, 0f), new PDist(0f, 0f), new PDist(0.25f, 0f), new PDist(5f, 0f), new PDist(0.5f, 0.5f),
                tunnelRadius, new PDist(1f, 0f), new PDist(1f, 0f));

        Collection<Structure> structures = structureDefinition.generateStructures(seed, region.getRegion());

        for (Structure structure : structures) {
            structure.generateStructure(new CaveStructureCallback(facet));
        }

        region.setRegionFacet(CaveFacet.class, facet);
    }
}
