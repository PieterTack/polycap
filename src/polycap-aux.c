/*
 * Copyright (C) 2018 Pieter Tack, Tom Schoonjans and Laszlo Vincze
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * */

#include "config.h"
#include "polycap-aux.h"
#include <stdlib.h>
#include <string.h>

char *polycap_strdup(const char *str) {
#ifdef HAVE__STRDUP
	return _strdup(str);
#elif defined(HAVE_STRDUP)
	return strdup(str);
#else
	char *dup= (char *)malloc( strlen(str)+1 );
	if (dup) strcpy(dup,str);
	return dup;
#endif
}

char *polycap_strndup(const char *str, size_t len) {
#ifndef HAVE_STRNDUP
	char *dup= (char *)malloc( len+1 );
	if (dup) {
		strncpy(dup,str,len);
		dup[len]= '\0';
	}
	return dup;
#else
	return strndup(str, len);
#endif
}
