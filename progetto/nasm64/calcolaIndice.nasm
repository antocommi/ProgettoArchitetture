section .data
section .bss
section .text

global calcolaIndice

;i equ 8	rdi
;j equ 12	rsi

calcolaIndice:
		push rbp
		mov	rbp, rsp
		;push rsi
		;push rdi	;inizio

		mov rax, rdi		;i
		dec rax				;i-1
		imul rax, rdi		;i*(i-1)
		sar rax, 1			;i*(i-1)/2
		add rax, rsi		;i*(i-1)/2+j

		;pop	rdi		;fine
		;pop	rsi
		mov	rsp, rbp
		pop	rbp
		ret