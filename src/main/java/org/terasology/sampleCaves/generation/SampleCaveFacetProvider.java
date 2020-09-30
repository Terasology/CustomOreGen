// Copyright 2020 The Terasology Foundation
// SPDX-License-Identifier: Apache-2.0
package org.terasology.sampleCaves.generation;

import org.terasology.customOreGen.PDist;
import org.terasology.customOreGen.Structure;
import org.terasology.customOreGen.StructureDefinition;
import org.terasology.customOreGen.VeinsStructureDefinition;
import org.terasology.entitySystem.Component;
import org.terasology.math.JomlUtil;
import org.terasology.nui.properties.Range;
import org.terasology.world.generation.ConfigurableFacetProvider;
import org.terasology.world.generation.FacetProvider;
import org.terasology.world.generation.GeneratingRegion;
import org.terasology.world.generation.Produces;

import java.util.Collection;

@Produces(SampleCaveFacet.class)
public class SampleCaveFacetProvider implements ConfigurableFacetProvider, FacetProvider {
    private long seed;
    private CaveFacetProviderConfiguration configuration = new CaveFacetProviderConfiguration();

    @Override
    public void setSeed(long seed) {
        this.seed = seed;
    }

    @Override
    public void process(GeneratingRegion region) {
        SampleCaveFacet facet = new SampleCaveFacet(region.getRegion(), region.getBorderForFacet(SampleCaveFacet.class));

        /* AnotherWorld
        PDist caveFrequency = new PDist(0.1f, 0f);
        PDist mainCaveRadius = new PDist(5f, 1f);
        PDist mainCaveYLevel = new PDist(32f, 32f);
        PDist tunnelLength = new PDist(45f, 10f);
        PDist tunnelRadius = new PDist(2f, 0.5f);
        */

        // note, try to contain the maxRange to an optimal range of neighboring chunks
        PDist caveFrequency = new PDist(configuration.frequency, configuration.frequency * 0.25f);
        PDist mainCaveRadius = new PDist(configuration.caveRadius, configuration.caveRadius * 0.5f);
        PDist mainCaveYLevel = new PDist(32f, 32f);
        PDist tunnelLength = new PDist(40f, 10f);
        PDist tunnelRadius = new PDist(configuration.tunnelRadius, configuration.tunnelRadius * 0.25f);

        StructureDefinition structureDefinition = new VeinsStructureDefinition(caveFrequency,
                mainCaveRadius, mainCaveYLevel, new PDist(4f, 1f), new PDist(0f, 0.1f), tunnelLength,
                new PDist(100f, 0f), new PDist(0f, 0f), new PDist(0.25f, 0f), new PDist(5f, 0f), new PDist(0.5f, 0.5f),
                tunnelRadius, new PDist(1f, 0f), new PDist(1f, 0f));

        Collection<Structure> structures = structureDefinition.generateStructures(seed, JomlUtil.from(region.getRegion()));

        for (Structure structure : structures) {
            structure.generateStructure(new SampleCaveStructureCallback(facet));
        }

        region.setRegionFacet(SampleCaveFacet.class, facet);
    }


    @Override
    public String getConfigurationName() {
        return "Caves";
    }

    @Override
    public Component getConfiguration() {
        return configuration;
    }

    @Override
    public void setConfiguration(Component configuration) {
        this.configuration = (CaveFacetProviderConfiguration) configuration;
    }

    private static class CaveFacetProviderConfiguration implements Component {
        @Range(min = 0, max = 1f, increment = 0.01f, precision = 2, description = "Cave Frequency")
        public float frequency = 0.1f;
        @Range(min = 0, max = 25f, increment = 1f, precision = 0, description = "Cave Radius")
        public float caveRadius = 8f;
        @Range(min = 0, max = 10f, increment = 1f, precision = 0, description = "Tunnel Radius")
        public float tunnelRadius = 4f;
    }
}
