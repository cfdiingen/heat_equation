        module variables
        integer i,j,k
        integer nx,ny,nz
        integer ptemp,fin,campo
        real alfax,alfay,alfaz
        real dx,dy,dz,dt
        real lx,ly,lz
        real fox,foy,foz,fmax
        real ti,ts,tn,te,toe,th,tb,psave
        real,allocatable,dimension(:,:,:)::T,T0
        real,allocatable,dimension(:)::x,y,z
        end module
        !------------------------------------------------------------
        program dif3d 
        use variables
        nx=30!i
        ny=30!i
        nz=30!i
        lx=1.0!f
        ly=1.0!f
        lz=1.0!f
        alfax=8.0E-6!f
        alfay=8.0E-6!f
        alfaz=8.0E-6!f
        dt=15.0!f
        fin=1500!i
        psave=4!i
        ti=100.0!f
        tn=20.0!f
        ts=20.0!f
        te=100.0!f
        toe=100.0!f
        th=100.0!f
        tb=100.0!f

        write(*,*)
        write(*,*)'------ Maha-dif -------'
        write(*,*)
        write(*,*)'------ Iniciando...-------'
        call alocar()
        call dif()
        call fourier()
        write(*,*)fmax


        if(fmax.le.0.25)then
        write(*,*)
        write(*,*)'------ Calculando...-------'
          call malla()
          call val_ini()
          call cf_y()
          call cf_x()
          call cf_z()
          call pasos()
          write(*,*)
          write(*,*)'------ Solucion :) ------'

        else
          write(*,*)
          write(*,*)'------ Error:Inestabilidad :| -------'
          write(*,*)
          write(*,*)'------ Todo Fourier debe ser < 0.25 -------'
        end if
          call info()
          call dealocar()
        end
        !------------------------------------------------------------
        subroutine alocar()
        use variables
        allocate(T(nx,ny,nz))
        allocate(T0(nx,ny,nz))
        allocate(x(nx))
        allocate(y(ny))
        allocate(z(nz))
        return
        end subroutine
        !------------------------------------------------------------
        subroutine dealocar()
        use variables
        deallocate(T)
        deallocate(T0)
        deallocate(x)
        deallocate(y)
        deallocate(z)
        return
        end subroutine
        !------------------------------------------------------------
        subroutine dif()
        use variables
        dx=lx/(nx-1)
        dy=ly/(ny-1)
        dz=lz/(nz-1)
        return
        end subroutine
        !------------------------------------------------------------
        subroutine fourier()
        use variables
        fox=dt*alfax/(dx**2.)
        foy=dt*alfay/(dy**2.)
        foz=dt*alfaz/(dz**2.)
        fmax=max(fox,foy,foz) 
        return
        end subroutine
        !------------------------------------------------------------
        subroutine malla()
        use variables
        do i=1,nx
          x(i)=dx*float(i-1) 
        end do
        do j=1,ny
          y(j)=dy*float(j-1) 
        end do
        do k=1,nz
          z(k)=dz*float(k-1) 
        end do
        return
        end subroutine
        !------------------------------------------------------------
        subroutine val_ini()
        use variables
        do i=1,nx
          do j=1,ny
            do k=1,nz
               T(i,j,k)=ti
            end do
          end do
        end do
        return
        end subroutine
        !------------------------------------------------------------
         
        subroutine cf_y()
        use variables
        do i=1,nx
          do k=1,nz
             T(i,1,k)=tn
             T(i,ny,k)=ts
           end do
        end do
        return
        end subroutine
        !------------------------------------------------------------
        subroutine cf_x()
        use variables
        do j=1,ny
          do k=1,nz
             T(1,j,k)=toe
             T(nx,j,k)=te
          end do
        end do
        return
        end subroutine
        !------------------------------------------------------------
        subroutine cf_z()
        use variables
        do i=1,nx
          do j=1,ny
             T(i,j,1)=tb
             T(i,j,nz)=th
          end do
        end do
        return
        end subroutine
        !------------------------------------------------------------
        subroutine pasos()
        use variables
        integer it
        it=1
        campo=1
        open(20,file="out3d.dat")

        do ptemp=1,fin
            call solucion()
            if(it.eq.psave)then
              call tech()
              write(*,*)'Guardando campo:',campo
              campo=campo+1
              it=1
            else
              it=it+1
            end if
        end do
        close(20)
        return
        end subroutine
        !------------------------------------------------------------
        subroutine solucion
        use variables
          do i=2,nx-1
             do j=2,ny-1
                do k=2,nz-1
                T0(i,j,k)=(T(i+1,j,k)-2.0*T(i,j,k)+T(i-1,j,k))*fox+&
                         (T(i,j+1,k)-2.0*T(i,j,k)+T(i,j-1,k))*foy+&
                         (T(i,j,k+1)-2.0*T(i,j,k)+T(i,j,k-1))*foz+&
                          T(i,j,k) 
                end do
             end do
          end do 
          do i=2,nx-1
             do j=2,ny-1
                do k=2,nz-1
                T(i,j,k)=T0(i,j,k)
                end do
             end do
          end do 
        return
        end subroutine
        !------------------------------------------------------------
       subroutine tech()
       use variables

       write(20,*)'variables="x", "y", "z", "T"'
       write(20,*)'zone I= ',nz,' J= ',ny,'K= ',nx,'datapacking=point'
       do i=1,nx
          do j=1,ny
             do k=1,nz
                write(20,*)x(i),y(j),z(k),T(i,j,k)
             end do
          end do
       end do
       return
       end subroutine
        !------------------------------------------------------------
       subroutine info
       use variables
       write(*,*)
       write(*,*)'Fourier(x) ------',fox
       write(*,*)'Fourier(y) ------',foy
       write(*,*)'Fourier(z) ------',foz
       write(*,*)'dx... -------',dx
       write(*,*)'dy... -------',dy
       write(*,*)'dz... -------',dz
       write(*,*)'dt... -------',dt
       write(*,*)
       write(*,*)'------ 2012 Viva Mexico Kbrns -|:)...-------'
       write(*,*)
       return
       end subroutine
        !------------------------------------------------------------
