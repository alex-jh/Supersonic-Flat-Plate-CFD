#pragma once
#include "State.h"
#include <memory>

namespace cavity_flow {
	class CavityFlow
	{
	public:
		CavityFlow();
		~CavityFlow();

		void Run();

	private:
		State state;
		Mesh mesh;
		std::vector<std::unique_ptr<Boundary> > boundaries;
	};
}

