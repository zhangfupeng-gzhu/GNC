



subroutine get_tnr_timescale_at_rh(tnr)
	use com_main_gw
	implicit none
	 real(8) tnr, rho, sma, lglambda,sigma, blkmass, alpha
	 real(8) nm2_tot
	 integer i
	 nm2_tot=0

	 do i=1, ctl%m_bins
		nm2_tot=nm2_tot+ctl%n0*ctl%asymptot(1,i)*ctl%bin_mass(i)**2
		!print*, "nm2tot, n0, asm, m=", nm2_tot, ctl%n0, ctl%asymptot(1,i), ctl%bin_mass(i)
	 end do
	 lglambda=log(mbh)
	 tnr=0.34*ctl%v0**3/nm2_tot/lglambda/2d0/pi/1d6
	 if(rid.eq.0)then
	 	print*,"TNR(rh)=",tnr, " Myr"
		!stop
	 end if
	 !stop
 end subroutine


subroutine output_sg_sample_track_txt(sp, fl)
	use com_main_gw
	implicit none
	character*(*) fl
	character*(100) str_flag
	integer i,n
	character*(20) tmpid
	type(particle_sample_type)::sp

	write(unit=tmpid,fmt="(I15)") sp%id
	open(unit=999,file=trim(adjustl(fl))//"_"//trim(adjustl(tmpid))//"_simple.txt")
	    print*, "fl=",trim(adjustl(fl))//"_simple.txt"
		write(unit=999, fmt="(A29, 4A20,A10, A30)") "time", "ac", "ec", "inc","om", "flag", "state"
        print*, "sp%length=",sp%length,sp%id
		do i=1, sp%length
            !print*, "sp%length=",sp%length
            !print*, "i=",i
			!call get_str_flag(sp%track(i)%state_flag, str_flag)
			!select case (sp%track(i)%state_flag)
			!case (state_inby_kl, state_inby_21_flyby,state_binary_MBH,&
			!		state_inby_21_exchange,state_bhb_mbh_inby_gw, &
			!		state_inby_gw_kl, state_inby_gw_3body,state_bhb_mbh_inby_gw_final,state_merge)
				write(unit=999, fmt="(E29.15, 4E20.10,I10)") sp%track(i)%time, sp%track(i)%ac, sp%track(i)%ec,&
            sp%track(i)%incout, sp%track(i)%omout, sp%track(i)%state_flag !, trim(adjustl(str_flag))
			!end select
			!print*, str_flag
		end do
	close(unit=999)
end subroutine

subroutine output_sams_sg_track_txt(sps, fl)
	use com_main_gw
	implicit none
	character*(*) fl
	integer i,n, num, num2
	type(particle_samples_arr_type)::sps
	!type(sample_type),pointer::sp
	character*(10) itmp, tmpid
	integer,parameter:: nummax=500
	character*(4) str

    print*, "track n=",sps%n
	num=0;num2=0
    if(.not.allocated(sps%sp))then
        print*, "sps_arr not allocated"
        return
    end if
	do i=1, sps%n
		associate(sp=>sps%sp(i))
		if(sp%write_down_track.ge.record_track_nes&
            .or.(ctl%trace_all_sample.ge.record_track_detail))then
            !print*, "i=",i, sp%exit_flag
			write(unit=itmp, fmt="(I10)") i+rid*nummax
			select case(sp%exit_flag)
			case(exit_boundary_max)
				num2=num2+1
				if(num2<nummax)then
					if(ctl%output_track_emri.ge.1)then
						select case(sp%obtype)
						case(star_type_ms)
							call output_sample_track_txt(sp,trim(adjustl(fl))//"/emri/sg/ms/evl_sg_"//trim(adjustl(itmp)))
						case(star_type_bh)

								call output_sample_track_txt(sp,trim(adjustl(fl))//"/emri/sg/bh/evl_sg_"//trim(adjustl(itmp)))
						
						end select
					end if
				end if
			case(exit_plunge_single)
				select case(sp%obtype)
				case(star_type_ms)
					!call output_sample_track_txt(sp,trim(adjustl(fl))//"/plunge/ms/evl_sg_"//trim(adjustl(itmp)))
				case(star_type_bh)
					if(ctl%output_track_plunge.ge.1)then
						call output_sample_track_txt(sp,trim(adjustl(fl))//"/plunge/bh/evl_sg_"//trim(adjustl(itmp)))
					end if
				end select
			case(exit_tidal_empty,exit_tidal_full)
				num=num+1
				if(num<nummax)then
					if(ctl%output_track_td.ge.1)then
						call output_sample_track_txt(sp,trim(adjustl(fl))//"/td/evl_sg_"//trim(adjustl(itmp)))
					end if
				end if
			case(exit_normal)
				!call output_sample_track_txt(sp,trim(adjustl(fl))//"/norm/evl_sg_"//trim(adjustl(itmp)))				
			end select
		end if
        end associate
	end do
end subroutine
