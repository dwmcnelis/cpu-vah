/*
 * HydroParameterTest.cpp
 *
 *  Created on: Aug 21, 2017
 *      Author: mjmcnelis
 */

#include <libconfig.h>
#include <unistd.h>

#include "edu/osu/rhic/harness/cli/catch.hpp"
#include "edu/osu/rhic/harness/hydro/HydroParameters.h"

TEST_CASE( "hydro parameters from file", "[harness][parameters]" ) {
	struct HydroParameters params;
	config_t config;
	config_init(&config);

	char *rootDirectory = NULL;
	size_t size;
	char pathToConfigFile[255];
	rootDirectory = getcwd(rootDirectory,size);
	sprintf(pathToConfigFile, "%s/rhic/rhic-harness/src/test/resources", rootDirectory);
	loadHydroParameters(&config, pathToConfigFile, &params);
	config_destroy(&config);

    REQUIRE( params.initialProperTimePoint == 0.5 );
    REQUIRE( params.shearViscosityToEntropyDensity == 0.2 );
}

TEST_CASE( "hydro default parameters", "[harness][parameters]" ) {
	struct HydroParameters params;
	config_t config;
	config_init(&config);
	loadHydroParameters(&config, "", &params);
	config_destroy(&config);

    REQUIRE( params.initialProperTimePoint == 0.1 );
    REQUIRE( params.shearViscosityToEntropyDensity == 0.0795775 );
}
