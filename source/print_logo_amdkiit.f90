subroutine print_logo_amdkiit()
USE version

write(6,"(A82)")      "       ____        ___    ___  _________    __    __  __  __  ____________ "
write(6,"(A82)")      "     /  __  \     |   \  /   ||   ____  \  |  |  /  /|  ||  ||____    ____|"
write(6,"(A82)")      "    /  /  \  \    |    \/    ||  |    \  \ |  | /  / |  ||  |     |  |     "
write(6,"(A82)")      "   /  /____\  \   |  |\  /|  ||  |     |  ||  |/  /  |  ||  |     |  |     "
write(6,"(A82)")      "  /  ________  \  |  | \/ |  ||  |     |  ||  |\  \  |  ||  |     |  |     "
write(6,"(A82)")      " /  /        \  \ |  |    |  ||  |____/  / |  | \  \ |  ||  |     |  |     "
write(6,"(A82)")      "/__/          \__\|__|    |__||_________/  |__|  \__\|__||__|     |__|     "
write(6,*)
!write(6,"(A)")repeat("*", 93)
write(6,"(A)")repeat("~", 93)
write(6, "(3A)")repeat(" ", 5),"A Plane Wave DFT Program for Ab Initio Dynamics at Higher Rungs of DFT Functionals",repeat(" ", 5)
write(6,"(A)")repeat("~", 93)
write(6,*)
write(6,"(4A)")repeat(" ", 35),"**  VERSION :",trim(project_ver),"  **"
write(6,*)
write(6,"(4A)")repeat(" ", 26),"Binary Build date :",compilation_time ,repeat(" ",38)
write(6,"(A)")repeat("*", 93)

end subroutine print_logo_amdkiit

