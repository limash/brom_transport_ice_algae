    program main
    
    !REVISION HISTORY:
    !Original author(s): Shamil Yakubov
    !dev_version
    
    use brom_transport, only: init_brom_transport, do_brom_transport, clear_brom_transport

    !initializing, writing from gotm data and fabm.yaml included
    call init_brom_transport()
    !main cycle
    call do_brom_transport()
    !clear all
    call clear_brom_transport()

    end program main
    
!-----------------------------------------------------------------------
! Copyright under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
