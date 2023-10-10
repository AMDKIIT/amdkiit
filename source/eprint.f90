subroutine eprint()
  USE kinds
  USE system_data_types
  write(6,"(A)")repeat("*", 93) 
  WRITE(6,'(1A,A30,A2,T38,F15.10,A5)')"#",'TOTAL ENERGY',"=",etotal,"A.U."
  WRITE(6,'(1A,A30,A2,T38,F15.10,A5)')"#",'KINETIC ENERGY',"=",eke,"A.U."
  WRITE(6,'(1A,A30,A2,T38,F15.10,A5)')"#",'ELECTROSTATIC ENERGY',"=",estat,"A.U."
  WRITE(6,'(1A,A30,A2,T38,F15.10,A5)')"#",'SELF ENERGY CONTRIBUTION',"=",dreal(eself),"A.U."
  WRITE(6,'(1A,A30,A2,T38,F15.10,A5)')"#",'OVERLAP ENERGY CONTRIBUTION',"=",esr,"A.U."
  WRITE(6,'(1A,A30,A2,T38,F15.10,A5)')"#",'LOCAL PP ENERGY',"=",dreal(eloc),"A.U."
  WRITE(6,'(1A,A30,A2,T38,F15.10,A5)')"#",'NON LOCAL PP ENERGY',"=",enl,"A.U."
  WRITE(6,'(1A,A30,A2,T38,F15.10,A5)')"#",'EXCHANGE CORRELATION ENERGY',"=",exc,"A.U."
  write(6,"(A)")repeat("*", 93)
end subroutine eprint

