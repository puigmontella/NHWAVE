










!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!
!    two_layer_slide.F
!
!    This file contains the subroutines which compute the slope-oriented, two-layer slide
!    based on granular flow rheology.  This does NOT contain Dmitri's viscous lower layer.
!
!    fluxes_ll
!    source_terms_ll
!    eval_huv_ll
!    convert_Ha_to_Hs
!
!    James Kirby, 6/27/16
!
!    Model described in:
!
!    Ma, G., Kirby, J. T., Hsu, T.-J. and Shi, F., 2015, "A two-layer granular landslide model for
!        tsunami wave generation: Theory and computation", Ocean Modelling, 93, 40-55, 
!        doi:10.1016/j.ocemod.2015.07.012 
!
!    This file is part of NHWAVE.
!
!    NHWAVE is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    NHWAVE is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with NHWAVE.  If not, see <http://www.gnu.org/licenses/>.!
!
!---------------------------------------------------------------------------------------
!
