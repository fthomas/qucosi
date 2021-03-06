// QuCoSi - Quantum Computer Simulation
// Copyright © 2009 Frank S. Thomas <f.thomas@gmx.de>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include <cppunit/ui/text/TestRunner.h>

#include <AlgorithmsTest.h>
#include <GateTest.h>
#include <QubitTest.h>
#include <VectorTest.h>

int main(int argc, char* argv[])
{
  CppUnit::TextUi::TestRunner runner;

  runner.addTest(QuCoSi::VectorTest::suite());
  runner.addTest(QuCoSi::QubitTest::suite());
  runner.addTest(QuCoSi::GateTest::suite());
  runner.addTest(QuCoSi::AlgorithmsTest::suite());
  runner.run();

  return 0;
}

// vim: shiftwidth=2 textwidth=78
