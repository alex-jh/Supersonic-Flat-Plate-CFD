#pragma once
#include "State.h"

namespace supersonic_nozzle {
	class SupersonicNozzle
	{
	public:
		SupersonicNozzle();
		~SupersonicNozzle();

		void Run();

	private:
		State state;
		Mesh mesh;
		std::vector<std::unique_ptr<Boundary> > boundaries;
	};
}

