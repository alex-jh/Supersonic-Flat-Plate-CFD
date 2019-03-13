#include <iostream>
#include <algorithm>
#include <sstream>
#include "SupersonicNozzle.h"
#include "CavityFlow.h"
using namespace std;


int main() {

	supersonic_nozzle::SupersonicNozzle problem;

	problem.Run();

	return 0;

}