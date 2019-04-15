        program iso
        
        implicit none
        integer*8 i,j,nmax,n(500),m(500),k,f
        integer*8 nn,mm,jj
        real*8 pf(500),ki(500),t(50)
        real*8 d(500,50)
        real*8 p(0:20,50)
        character*20 pfact,kint,ass,times

        read(*,*)f
        read(*,*)ass
        read(*,*)kint
        read(*,*)pfact
        read(*,*)times
        open(89,file=ass)
        open(90,file=pfact)
        open(91,file=kint)
        open(92,file=times)
        open(99,file="out.Dpred")

        nmax=0
        do j=1,400
          ki(j)=0.0
          pf(j)=0.0
        end do
        do i=1,50
          t(i)=-1.0
        end do
        do j=1,400
          read(89,*,END=90)i,nn,mm
          n(j)=nn
          m(j)=mm
        end do
90      continue
        do j=1,400
          read(90,*,END=100)i,pf(i)
          read(91,*,END=100)i,ki(i)
          if (ki(i) .gt. 0.0) nmax=nmax+1
          write(*,*)i,pf(i),ki(i)
        end do
100     continue
        write(*,*)f,n(f),m(f),nmax
        do i=n(f),m(f)
          write(*,*)i,pf(i),ki(i)
        end do
        do i=1,50
          read(92,*,END=110)t(i)
        end do
110        continue
        do i=1,50
           if (t(i).gt.0.0) then
              jj=0
              do j=n(f),m(f)
                 if (ki(j).gt.0.0) then
                    jj=jj+1
                    d(jj,i)=1.0-exp(-ki(j)*60*t(i)/exp(pf(j)))
                 end if
              end do
           end if
        end do
        write(*,*)n(f),m(f),jj
        do i=1,50
           if (t(i).gt.0.0) then
              write(99,'(99(g12.5,1x))')t(i),(d(j,i),j=1,jj)
           end if
        end do        

        end
