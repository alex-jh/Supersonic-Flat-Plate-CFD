#include "supersonic_flat_plate.h"
#include "supersonic_cone.h"
#include "supersonic_plate.h"
#include "supersonic_rocket_nozzle.h"

void RunSupersonicFlatPlate() {
    using namespace supersonic_flat_plate;
    SupersonicFlatPlate supersonic_flat_plate;
    supersonic_flat_plate.Run();
}

void RunSupersonicCone() {
    using namespace supersonic_cone;
    SupersonicCone supersonic_cone;
    supersonic_cone.Run();
}

void RunSupersonicPlate() {
	using namespace supersonic_plate;
	SupersonicPlate supersonic_plate;
	supersonic_plate.Run();
}

void RunSupersonicRocketNozzle() {
	using namespace supersonic_rocket_nozzle;
	SupersonicRocketNozzle supersonic_rocket_nozzle;
	supersonic_rocket_nozzle.Run();
}

int main() {
	RunSupersonicRocketNozzle();

	return 0;
}