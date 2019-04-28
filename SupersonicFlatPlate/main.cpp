#include "supersonic_flat_plate.h"
#include "supersonic_cone.h"
#include "supersonic_plate.h"

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

int main() {
	RunSupersonicCone();

	return 0;
}