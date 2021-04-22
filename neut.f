************************************************************************
*      Programa para estudiar el comportamiento de un haz de neutrones
*      al interaccionar con material de tipo org nico.
************************************************************************
       include 'ran2_r8.f'
       implicit real*8 (a-h,m,l,o-z)
       parameter (neu=10000)
       parameter (j=1200)
       dimension prof(j),cont(j)

       pi=4.0d0*atan(1.0d0)
       iseed=-23456777

       !masas y secciones eficaces  (mp es masa del proyectil empleado, no prot¢n)
       mp=1.0086649156d0; mH=1.00797d0; mC=12.0107d0
       mO=15.9994d0; mN=14.0067d0; uH=0.777d0
       uC=0.041d0; uO=0.1d0; uN=0.00555d0
       u=uH+uC+uO+uN
       
       !apartado II, profundidad
       pmin=0.0d0  !si ponemos 0 el apartado I
       pmax=25.0d0
       dp=(pmax-pmin)/dble(j)

       do i=1,j
           prof(i)=pmin+i*dp
           cont(i)=0
       end do
       
       open(10,file='xyz.dat')
       open(20,file='prof1000.dat')
       do n=1,neu
         E0=1000.0d0 !E_0 (KeV)
         Emin=1.0d-1 !energ¡a de termalizaci¢n

         !Condiciones iniciales
         x1=0.0d0; x2=0.0d0; x3=0.0d0; O=pi/2.0d0; Y=pi/2.0d0
         write(10,*) n,[x1, x2, x3]
         L=-log(ran2(iseed))/u !distancia sin interacci¢n
         x3=L
         write(10,*) n,[x1, x2, x3]  !x,y,z
         do while (E0.gt.Emin)
           L=-log(ran2(iseed))/u
           R=ran2(iseed) !elecci¢n del choque (de mayor prob. a menor)
           if (R.le.0.842d0) then
             mr=mH; go to 100
           else if (R.gt.0.842d0.and.R.le.0.950d0) then
             mr=mO; go to 100
           else if (R.gt.0.950d0.and.R.le.0.994d0) then
             mr=mC; go to 100
           else if (R.gt.0.994d0) then
             mr=mN
           end if

           !Aqu¡ comienza el c lculo de la dispersi¢n
         
  100      xi=acos(1.0d0-2.0d0*ran2(iseed))  ! ngulos aleatorios
           eta=2.0d0*pi*ran2(iseed)

           ! ngulos SCM
           xirec=pi-xi
           etarec=eta+pi
           
           ! ngulos SR del neutr¢n
           the=atan(1.0d0/(mp/(mr*sin(xi))+1.0d0/tan(xi)))
           phi=eta
           therec=xirec/2.0d0
           phirec=etarec

           !energ¡a perdida en el choque
           Er=E0*4.0d0*mp*mr*(cos(therec)/(mp+mr))**2
           E0=E0-Er

           ! ngulos SRL  (expresi¢n compacta)
           a=sin(the);b=cos(the);c=sin(phi);d=cos(phi);e=cos(O);f=sin(O)
           cosOrec=-a*d*f+b*e
           Orec=acos(cosOrec)
           cosYrec=(a*d*e*cos(Y)-a*c*sin(Y)+b*f*cos(Y))/sin(Orec)
           O=Orec
           Y=acos(cosYrec)

           !nuevas coordenadas x,y,z tras el choque
           x1=x1+L*cos(O)
           x2=x2+L*sin(O)*cos(Y)
           x3=x3+L*sin(O)*sin(Y)

           write (10,*) n,[x1, x2, x3] !x,y,z
         end do
         write(10,*)  !l¡nea en blanco para separar neutrones

*        apartado II, profundidad alcanzada
         do i=1,j-1
           if(x3.ge.prof(i).and.x3.lt.prof(i+1)) then
            cont(i)=cont(i)+1
           end if
         end do
       end do

       do i=1,j !bucle para almacenar profundidades
         write(20,*) prof(i),dble(cont(i))
       end do

       close(10)
       close(20)
       write(*,*) 'Programa finalizado correctamente. Pulse una tecla.'
       stop
       end
