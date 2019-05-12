c---  dimension 8 operators
      Adim8(1,2,2)= A1pp(za,zb)*im*(four,zero)
      Adim8(1,2,1)= A1pm(za,zb)*im*ctwo
      Adim8(1,1,2)= conjg(Adim8(1,2,1))
      Adim8(1,1,1)= conjg(Adim8(1,2,2))

      Adim8(2,2,2)= A2pp(za,zb)*im*(8._dp,zero)
      Adim8(2,2,1)= czip
      Adim8(2,1,2)= czip
      Adim8(2,1,1)= conjg(Adim8(2,2,2))

      Adim8(3,2,2)= A3pp(za,zb)*im*(8._dp,zero)
      Adim8(3,2,1)= A3pm(za,zb)*im*(8._dp,zero)
      Adim8(3,1,2)= conjg(Adim8(3,2,1))
      Adim8(3,1,1)= conjg(Adim8(3,2,2))

      Adim8(4,2,2)= A4pp(za,zb)*im*(four,zero)
      Adim8(4,2,1)= A4pm(za,zb)*im*(four,zero)
      Adim8(4,1,2)= conjg(Adim8(4,2,1))
      Adim8(4,1,1)= conjg(Adim8(4,2,2))

      Adim8(5,2,2)= A5pp(za,zb)*im*wmass**2
      Adim8(5,2,1)= A5pm(za,zb)*im*wmass**2
      Adim8(5,1,2)= conjg(Adim8(5,2,1))
      Adim8(5,1,1)= conjg(Adim8(5,2,2))

      Adim8(6,2,2)= A6pp(za,zb)*im*(four,zero)*wmass**2
      Adim8(6,2,1)= czip
      Adim8(6,1,2)= czip
      Adim8(6,1,1)= conjg(Adim8(6,2,2))

      do h1=1,2
         do h2=1,2
            Adim8(:,h1,h2)=Adim8(:,h1,h2)*cdim8(:)
         enddo
      enddo