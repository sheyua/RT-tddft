This repository contains the Real Time-Time Dependent Density Functional
Theory (rt-tddft) extension which enables
and the Quantum Espresso code on top of it.

To compile rt-tddft, please compile the quantum espresso code (QE-5.2.0)
first.

	cd rt-tddft
	autoconf
	./configure --with-qepath=[YOUR QE PREFIX]
	make build

#* RT-tddft README

This repository simulate how real-time electronic states evolve after unsetting
and initial bias voltage between left and right side of a nano-capacity like / 
nano-junction like structure.

## Ground State Calculation

	Biased ground state are computed by the pwscf routine, fictious dipole
	moment is cancelled with a mirror image of the computational cell. please
	see the graphene example for how to setup the bias voltage

## Time Propagation

	The initial bias voltage can be unset with a linear decay form. Time 
	integration is carried out with the following schemes which all strictly
	conserve the charge:

		'CN'    : Crank-Nicolson, O(dt^2) local error, unconditionally stable
		'CN2'   : second order Crank-Nicolson, O(dt^3) local error
		'CN-mid': mid point Crank-Nicolson, O(dt^3) local error, unconditionally
		stable

	The integration schemes are driven by two different linear solvers:

		'itsolver': Iterative solver
		'cgsolver': Conjugate Gradient square solver, implemented by Xiaofeng
		Qian at MIT

#* Quantum Espresso README

This is the distribution of the Quantum ESPRESSO suite of codes (ESPRESSO: 
opEn-Source Package for Research in Electronic Structure, Simulation, 
and Optimization), promoted by the IOM-DEMOCRITOS National Simulation Center 
of the Italian CNR (http://www.democritos.it). 

Quick installation instructions for the impatient:

	unzip espresso.zip
	cd espresso
	./configure [options]
	make pw

("make" alone prints a list of acceptable targets). Binaries go in bin/.
For more information, see the general documentation in directory Doc/, 
package-specific documentation in */Doc/, and the web site
http://www.quantum-espresso.org/

All the material included in this distribution is free software;
you can redistribute it and/or modify it under the terms of the GNU
General Public License as published by the Free Software Foundation;
either version 2 of the License, or (at your option) any later version.

These programs are distributed in the hope that they will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
675 Mass Ave, Cambridge, MA 02139, USA.
