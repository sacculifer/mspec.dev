!!------------------------------------------------------------------------------
!! mspec_zoo.F90 is part of mspec model developed at KSE.
!! First version (21.05.2019) by Ovidio Garcia ovidio.garcia@hzg.de
!! This file includes code related with zooplankton modelling using trait approach.
!! Also, it includes simple_zoo, a parametrization of zooplankton module.
!!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
subroutine grazing_forcing(self,Phy,F_T,Mean,zoo_pref,cop_pref,Zoo,grazing,Lz,I_max)
  !! grazing forcing implemented by Taherdazah et al. 2019.

  implicit none
  class (type_hzg_mspec),intent(in) :: self
  real(rk),dimension(self%phyto_num),intent(in) :: Phy ! Phytoplankton biomass concentration, mmol-C m^-3
  real(rk),dimension(self%zoo_num),intent(in) :: Zoo, Lz !Zooplankton biomass
  real(rk),dimension(2),intent(in) :: F_T  !(2) Temperature dependency for zooplankton
  real(rk),dimension(self%zoo_num,self%phyto_num),intent(in) :: zoo_pref
  real(rk),dimension(self%num_ciliat), intent(in) :: cop_pref
  real(rk),intent(in) :: Mean  ! Community mean cell size, log_e ESD (mu m)
  real(rk),dimension(self%zoo_num), intent(out) :: I_max
  real(rk),dimension(self%phyto_num),intent(out) :: grazing !Grazing forcing, mmol-C m^-3 d^-1

  real(rk),dimension(self%zoo_num,self%phyto_num) :: eff_food,mean_eff_food,graz_j,graz_j_spec
  real(rk),dimension(self%zoo_num) :: eff_food_con,mean_eff_food_con,glob_graz
  real(rk) :: eff_food_ciliat,nom,dom,a_zoo,I_max_star, x,g_x,G_min,mean_eff_food_ciliat,G,downreg
  integer :: i,j,c

  grazing = 0.0_rk
  eff_food_ciliat = 0.0_rk
  eff_food = 0.0_rk
  mean_eff_food_ciliat = 0.0_rk
  mean_eff_food = 0.0_rk
  I_max = 0.0_rk
  G_min = self%a_gr*self%R_A/self%y
  graz_j = 0.0_rk
  glob_graz = 0.0_rk

  !! Effective food concentration
  do c=1,self%num_ciliat
    eff_food_ciliat=eff_food_ciliat+(Zoo(c)*cop_pref(c))
  end do

  do j=1,self%zoo_num
    do i=1,self%phyto_num
      eff_food(j,i)=Phy(i)*zoo_pref(j,i)
    end do
    eff_food_con(j)=sum(eff_food(j,:))
  end do

  !! effective food con for copepoda-> Copepoda graze on phytoplankton + ciliates
  eff_food_con(self%zoo_num) = eff_food_con(self%zoo_num)+eff_food_ciliat

  !! Average effective food concentration
  do c=1,self%num_ciliat
    mean_eff_food_ciliat=mean_eff_food_ciliat+Zoo(c)*cop_pref(c) * Lz(c)
  end do

  do j=1,self%zoo_num
    do i=1,self%phyto_num
      mean_eff_food(j,i)=Phy(i)*zoo_pref(j,i)*self%log_ESD(i)
    end do
    mean_eff_food_con(j)=sum(mean_eff_food(j,:))
  end do

  mean_eff_food_con(self%zoo_num)=mean_eff_food_con(self%zoo_num)+mean_eff_food_ciliat

  !! Maximum ingestion rate
  do j=1,self%zoo_num
    a_zoo=self%a_Im0*(Lz(j)+self%Lz_star(j))
    I_max_star=self%I_max0*F_T(2)*exp(a_zoo+(2.0_rk-a_zoo)*self%Lz_star(j)+(a_zoo-3.0_rk)*Lz(j))
    I_max(j)=I_max_star*exp(-self%sel(j)*(self%Lz_star(j)-mean_eff_food_con(j)/eff_food_con(j))**2)
  end do

  !! Calculation of grazing, x: food processing ratio, g_x: functional response
  do j=1,self%zoo_num
    x=self%a_gr * eff_food_con(j)/I_max(j)
    g_x=1.0_rk-exp(-x) !((1.0_rk-x**self%n_syn)/(1.0_rk-x**(self%n_syn+1.0_rk)))*x
    glob_graz(j)=self%graz_const*Zoo(j)*I_max(j)*g_x
  end do

  !! Threshold
  do i=1,self%phyto_num
    do j=1,self%zoo_num
      !graz_j_spec(j,i)=zoo_pref(j,i)*Phy(i)/(eff_food_con(j)+eps)
      graz_j(j,i)=glob_graz(j)*zoo_pref(j,i)*Phy(i)/(eff_food_con(j)+eps)
    end do
  end do

  do j=1,self%zoo_num
    G = sum(graz_j(j,:))
    downreg=1.0_rk/(1.0_rk+exp(-(G-G_min)/0.05_rk))
    graz_j(j,:)=downreg*graz_j(j,:)
  end do

  do i=1,self%phyto_num
    grazing(i)=sum(graz_j(:,i))
  end do

  return
end subroutine grazing_forcing

!!------------------------------------------------------------------------------
subroutine dummy_grazing(self,Phy,grazing)
  !! Calculates phytoplankton losses due grazing using a squared loss rate.
  !! Intent to be used when zooplankton grazing is turned off.
  !! OG 21.05.2019

  implicit none
  class (type_hzg_mspec),intent(in) :: self
  real(rk),dimension(self%phyto_num),intent(in) :: Phy ! Phytoplankton biomass concentration, mmol-C m^-3
  real(rk),dimension(self%phyto_num),intent(out) :: grazing ! Grazing forcing, mmol-C m^-3 d^-1

  grazing = 0.04_rk*Phy**2

  return
end subroutine dummy_grazing

!!------------------------------------------------------------------------------
