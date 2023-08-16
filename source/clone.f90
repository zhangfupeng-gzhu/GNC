
logical function Istransit(en0,en1,i)
	use com_main_gw
    implicit none
	real(8) en0, en1, en0p, en1p
	integer i,jk,leveln
	i=-1
	Istransit=.false.
!	level0=log10(en0);level1=log10(en1)
    
	leveln=log10clone_emax
	!print*, leveln
	!stop
    en0p=log10(en0/ctl%clone_e0)
    en1p=log10(en1/ctl%clone_e0)
    !if(ctl%chattery.ge.3) then
    !    print*, "clone_emax=",clone_emax
    !    print*, "en0p, en1p=", en0p, en1p, leveln
    !end if
	do jk=1, leveln
    !    if(ctl%chattery.ge.3) print*, "jk=",jk 
		if(en0p<jk.and.en1p>jk)then
			Istransit=.true.
			i=jk
			return 
		end if
	end do
    !if(ctl%chattery.ge.3)then
    !    print*, "istransit=",istransit
    !endif
!	print*, "error in istransit",en0,en1
!	stop
end function

logical function Invtransit(en0,en1,i)
	use com_main_gw
    implicit none
	real(8) en0, en1, en0p, en1p
	integer i,jk,leveln

	Invtransit=.false.
	i=-1

	leveln=log10clone_emax
    en0p=log10(en0/ctl%clone_e0)
    en1p=log10(en1/ctl%clone_e0)
    if(ctl%chattery.ge.4)then
        print*, "en0p, en1p, clone_e0, clone_emax=", en0p, en1p, ctl%clone_e0, clone_emax
    end if
	do jk=1, leveln
        !if(ctl%chattery.ge.3)then
        !    print*, "jk=", jk, en0p>10**jk, en1p<10**jk
        !end if
		if(en0p>jk.and.en1p<jk)then
			Invtransit=.true.
			i=jk
			return 
		end if
	end do
!	print*, "error in invtransit",en0,en1
!	stop
end function

subroutine create_init_clone_particle(pt,spen0,time)
	use com_main_gw
	implicit none
	type(chain_pointer_type)::pt
	real(8) e0,spen0,time
	integer nlvl, i
    
	nlvl=int(log10(spen0/ctl%clone_e0))
	if(nlvl.ge.1)then
        !if(ctl%chattery.ge.3)then
        !    print*, "clone particle generated"
        !end if
        select type (ca=>pt%ob)
        type is(particle_sample_type)
            call get_mass_idx(ca%m, sample_mass_idx)
		    call create_clone_particle(pt,nlvl,ctl%clone_factor(sample_mass_idx)**nlvl, time)
        !type is(sample_type)
		!	call get_mass_idx(ca%m, sample_mass_idx)
        !    call create_clone_particle(pt,nlvl,ctl%clone_factor(sample_mass_idx)**nlvl, time)
        end select
	end if
end subroutine

subroutine clone_scheme(pt, en0, en1, amplifier,time, out_flag_clone)
    use com_main_gw
    implicit none
    type(chain_pointer_type)::pt
    integer amplifier
    real(8) en0, en1, time
    integer lvl, out_flag_clone
    logical Istransit, Invtransit
    real(8) rnd, tmp
    out_flag_clone=0
    if(Istransit(en0, en1, lvl))then
        if(ctl%chattery.ge.4)then
            if(ctl%ntasks>1)then
                write(unit=chattery_out_unit,fmt=*) "time,en0,en1=",time/1d6/(2*pi), en0,en1
                write(unit=chattery_out_unit,fmt=*) "crossing ",lvl,", create clone particle",  pt%idx
            else
                print*, "time,en0,en1=",time/1d6/(2*pi), en0,en1
                print*, "crossing ",lvl,", create clone particle",  pt%idx
            end if
        end if
        !print*, "ac_b:pt%ed%idx=",pt%ed%idx
        call create_clone_particle(pt,lvl,amplifier,time)
        !print*, "ac_e:pt%ed%idx=",pt%ed%idx
        !read(*,*)
        !Iscloned(lvl)=.true.
        out_flag_clone=1
        if(ctl%chattery.ge.4) write(*,*) "-------clone particle--",lvl,pt%idx, pt%ed%idx
    end if	
    !if(ctl%chattery.ge.3) print*, "???????", ctl%del_cross_clone
    if(ctl%del_cross_clone.ge.1)then
        !if(ctl%chattery.ge.3)then
        !    print*," nhiar, lvl=", sample%source,lvl
        !end if
        if(Invtransit(en0, en1,lvl))then
            tmp=rnd(0d0,dble(amplifier))
            if(tmp>1d0)then
                if(ctl%chattery.ge.4)then
                    if(ctl%ntasks>1)then
                        write(unit=chattery_out_unit,fmt=*) "en0,en1=",en0,en1
                        write(unit=chattery_out_unit,fmt=*) &
                            "Invcrossing, deleting clone particle",  pt%idx
                    else
                        print*, "en0,en1=",en0,en1
                        print*, "Invcrossing, deleting clone particle",  pt%idx
                    end if
                end if
            !	read(*,*)
                out_flag_clone=100
                return
            end if
        end if		
    end if
end subroutine
