/*
 * InitialConditionParameterTest.cpp
 *
 *  Created on: Oct 29, 2015
 *      Author: bazow
 */

#include <libconfig.h>
#include <unistd.h>

#include "edu/osu/rhic/harness/cli/catch.hpp"
#include "edu/osu/rhic/harness/ic/InitialConditionParameters.h"

TEST_CASE( "ic parameters from file", "[harness][parameters]" ) {
	struct InitialConditionParameters params;
	config_t config;
	config_init(&config);

	char *rootDirectory = NULL;
	size_t size;
	char pathToConfigFile[255];
	rootDirectory = getcwd(rootDirectory,size);
	sprintf(pathToConfigFile, "%s/rhic/rhic-harness/src/test/resources", rootDirectory);
	loadInitialConditionParameters(&config, pathToConfigFile, &params);
	config_destroy(&config);

    REQUIRE( params.initialConditionType == 0 );
    REQUIRE( params.numberOfNucleonsPerNuclei == 63 );
    REQUIRE( params.initialEnergyDensity == 0.6 );
    REQUIRE( params.scatteringCrossSectionNN == 52.0 );
    REQUIRE( params.impactParameter == 20.0 );
    REQUIRE( params.fractionOfBinaryCollisions == 0.1 );
    REQUIRE( params.rapidityVariance == 0.2 );
    REQUIRE( params.rapidityMean == 0.3 );
}

TEST_CASE( "ic default parameters", "[harness][parameters]" ) {
	struct InitialConditionParameters params;
	config_t config;
	config_init(&config);
	loadInitialConditionParameters(&config, "", &params);
	config_destroy(&config);

    REQUIRE( params.initialConditionType == 2 );
    REQUIRE( params.numberOfNucleonsPerNuclei == 208 );
    REQUIRE( params.initialEnergyDensity == 1 );
    REQUIRE( params.scatteringCrossSectionNN == 62 );
    REQUIRE( params.impactParameter == 7 );
    REQUIRE( params.fractionOfBinaryCollisions == 0.5 );
    REQUIRE( params.rapidityVariance == 0.5 );
    REQUIRE( params.rapidityMean == 0.5 );
}
