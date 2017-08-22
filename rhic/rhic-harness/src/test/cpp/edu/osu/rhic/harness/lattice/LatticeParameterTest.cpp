/*
 * LatticeParameterTest.cpp
 *
 *  Created on: Oct 29, 2015
 *      Author: bazow
 */

#include <libconfig.h>
#include <unistd.h>

#include "edu/osu/rhic/harness/cli/catch.hpp"
#include "edu/osu/rhic/harness/lattice/LatticeParameters.h"

TEST_CASE( "lattice parameters from file", "[harness][parameters]" ) {
	struct LatticeParameters params;
	config_t config;
	config_init(&config);

	char *rootDirectory = NULL;
	size_t size;
	char pathToConfigFile[255];
	rootDirectory = getcwd(rootDirectory,size);
	sprintf(pathToConfigFile, "%s/rhic/rhic-harness/src/test/resources", rootDirectory);
	loadLatticeParameters(&config, pathToConfigFile, &params);
	config_destroy(&config);

    REQUIRE( params.numLatticePointsX == 10 );
    REQUIRE( params.numLatticePointsY == 15 );
    REQUIRE( params.numLatticePointsRapidity == 20 );
    REQUIRE( params.numProperTimePoints == 133 );
    REQUIRE( params.latticeSpacingX == 0.1 );
    REQUIRE( params.latticeSpacingY == 0.2 );
    REQUIRE( params.latticeSpacingRapidity == 1 );
    REQUIRE( params.latticeSpacingProperTime == 0.5 );
}

TEST_CASE( "lattice default parameters", "[harness][parameters]" ) {
	struct LatticeParameters params;
	config_t config;
	config_init(&config);
	loadLatticeParameters(&config, "", &params);
	config_destroy(&config);

    REQUIRE( params.numLatticePointsX == 128 );
    REQUIRE( params.numLatticePointsY == 128 );
    REQUIRE( params.numLatticePointsRapidity == 64 );
    REQUIRE( params.numProperTimePoints == 10 );
    REQUIRE( params.latticeSpacingX == 0.08 );
    REQUIRE( params.latticeSpacingY == 0.08 );
    REQUIRE( params.latticeSpacingRapidity == 0.3 );
    REQUIRE( params.latticeSpacingProperTime == 0.01 );
}

